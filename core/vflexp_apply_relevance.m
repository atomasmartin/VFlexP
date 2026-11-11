function out = vflexp_apply_relevance(session, threshold, stateRel, stateMeta, dryRun)
%VFLEXP_APPLY_RELEVANCE Plan or apply model-complexity reductions from state relevance.
%
%   applied = vflexp_apply_relevance(session, threshold)
%   planTable = vflexp_apply_relevance(session, threshold, stateRel, stateMeta, true)
%
%   Simplification rules (when ALL states of an algebraic-capable sub-block
%   satisfy relevance < threshold):
%     - AC line:  dynamic -> static omega
%     - AC load:  series RL dynamic -> series RL static omega
%                 parallel RL dynamic -> parallel RL static omega
%                 dynamic PQ -> PQ
%     - DC line:  DC dynamic -> DC static
%     - DC load:  DC DynamicP -> DC ConstP
%     - MGES/SG:  order 8 -> order 4 (when stator_emt states < threshold;
%                 if rotor_flux states also < threshold -> order 2)
%     - VSM/GFM/intconv_VSM/intconv_GFM: Dynamic VCC -> Static VCC
%     - GFL/intconv_GFL: Dynamic CC -> Static CC; Dynamic PLL -> Static PLL
%     - AC/DC buses: dynamic -> static omega / DC static (evaluated last,
%                 guarded against algebraic-loop blockers)

if nargin < 5 || isempty(dryRun)
    dryRun = false;
end

if nargin < 2 || isempty(threshold)
    if isfield(session.config, 'relevance') && isfield(session.config.relevance, 'threshold')
        threshold = double(session.config.relevance.threshold);
    else
        threshold = 1e-3;
    end
end

if nargin < 3 || isempty(stateRel) || nargin < 4 || isempty(stateMeta)
    if isempty(session.relevanceResult)
        vflexp_relevance(session);
    end
    stateRel  = session.relevanceResult.stateRelevance;
    stateMeta = session.relevanceResult.stateMeta;
end

stateRel = double(stateRel(:));
nStates  = numel(stateRel);

if nStates == 0
    planTbl = local_empty_plan_table();
    if dryRun
        out = planTbl;
    else
        out = struct('nSimplified', 0, 'threshold', threshold, 'planTable', planTbl, 'blocked', strings(0,1));
    end
    return;
end

below = stateRel < threshold;
planRows = cell(0, 8);
blockedMsgs = strings(0, 1);

% Track unit-level complexity changes cleanly (Note: Simulink unit masks apply
% popup complexity across instances of that unit type, so we require ALL
% instances of a unit type to satisfy < threshold before simplifying its popup)
pendingUnitOverrides = struct();

% 1. AC Lines (per-branch)
acLineIdx = find(stateMeta.ElementKind == "ac_line" & stateMeta.SubGroup == "network");
acLineRows = unique(stateMeta.TableRow(acLineIdx), 'stable');
for r = acLineRows(:)'
    if r <= 0 || r > height(session.tables.Line), continue; end
    m = acLineIdx(stateMeta.TableRow(acLineIdx) == r);
    if all(below(m))
        oldC = char(strtrim(session.tables.Line.Complexity(r)));
        newC = 'static omega';
        if ~strcmpi(oldC, newC) && ~strcmpi(oldC, 'static')
            planRows(end+1, :) = {"ac_line", char(stateMeta.ElementId(m(1))), ...
                string(oldC), string(newC), max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
        end
    end
end

% 2. DC Lines (per-branch)
dcLineIdx = find(stateMeta.ElementKind == "dc_line" & stateMeta.SubGroup == "network");
dcLineRows = unique(stateMeta.TableRow(dcLineIdx), 'stable');
for r = dcLineRows(:)'
    if r <= 0 || r > height(session.tables.Line), continue; end
    m = dcLineIdx(stateMeta.TableRow(dcLineIdx) == r);
    if all(below(m))
        oldC = char(strtrim(session.tables.Line.Complexity(r)));
        newC = 'DC static';
        if ~strcmpi(oldC, newC)
            planRows(end+1, :) = {"dc_line", char(stateMeta.ElementId(m(1))), ...
                string(oldC), string(newC), max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
        end
    end
end

% 3. AC Loads (per-bus)
acLoadIdx = find(stateMeta.ElementKind == "ac_load" & stateMeta.SubGroup == "network");
acLoadRows = unique(stateMeta.TableRow(acLoadIdx), 'stable');
for r = acLoadRows(:)'
    if r <= 0 || r > height(session.tables.Bus), continue; end
    m = acLoadIdx(stateMeta.TableRow(acLoadIdx) == r);
    if all(below(m))
        oldC = char(strtrim(session.tables.Bus.Complexity_load(r)));
        oldLow = lower(oldC);
        if contains(oldLow, 'parallel rl')
            newC = 'parallel RL static omega';
        elseif contains(oldLow, 'dynamic pq')
            newC = 'PQ';
        else
            newC = 'series RL static omega';
        end
        if ~strcmpi(oldC, newC)
            planRows(end+1, :) = {"ac_load", char(stateMeta.ElementId(m(1))), ...
                string(oldC), string(newC), max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
        end
    end
end

% 4. DC Loads (per-bus)
dcLoadIdx = find(stateMeta.ElementKind == "dc_load" & stateMeta.SubGroup == "network");
dcLoadRows = unique(stateMeta.TableRow(dcLoadIdx), 'stable');
for r = dcLoadRows(:)'
    if r <= 0 || r > height(session.tables.Bus), continue; end
    m = dcLoadIdx(stateMeta.TableRow(dcLoadIdx) == r);
    if all(below(m))
        oldC = char(strtrim(session.tables.Bus.Complexity_load(r)));
        newC = 'DC ConstP';
        if ~strcmpi(oldC, newC)
            planRows(end+1, :) = {"dc_load", char(stateMeta.ElementId(m(1))), ...
                string(oldC), string(newC), max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
        end
    end
end

% 5. Synchronous Machines (MGES, SG): order 8 (EMT8) -> order 4 (RMS4)
%    Rotor fluxes (psifd, psikd1, psikq1, psikq2) have no algebraic form and are kept.
for uCell = {'MGES', 'SG'}
    uType = uCell{1};
    stIdx = find(stateMeta.UnitType == uType & stateMeta.SubGroup == "stator_emt");
    if ~isempty(stIdx) && all(below(stIdx))
        curStator = local_get_unit_mode(session, uType, 'dynamic_stator', 'order 8');
        if strcmpi(curStator, 'order 8') || strcmpi(curStator, 'order 8 nomatrix')
            newStator = 'order 4';
            pendingUnitOverrides.(uType).dynamic_stator = newStator;
            planRows(end+1, :) = {"machine", uType, string(curStator), string(newStator), ...
                max(stateRel(stIdx)), strjoin(stateMeta.Label(stIdx), ', '), 0, "dynamic_stator"}; %#ok<AGROW>
        end
    end
end

% 6. Grid-Forming Converters (VSM, GFM, intconv_VSM, intconv_GFM): VCC Dynamic -> Static
for uCell = {'VSM', 'GFM', 'intconv_VSM', 'intconv_GFM'}
    uType = uCell{1};
    vccIdx = find(stateMeta.UnitType == uType & stateMeta.SubGroup == "vcc");
    if ~isempty(vccIdx) && all(below(vccIdx))
        curVcc = local_get_unit_mode(session, uType, 'dynamic_VCC', 'Dynamic');
        if ~strcmpi(curVcc, 'Static')
            pendingUnitOverrides.(uType).dynamic_VCC = 'Static';
            planRows(end+1, :) = {"converter_vcc", uType, string(curVcc), "Static", ...
                max(stateRel(vccIdx)), strjoin(stateMeta.Label(vccIdx), ', '), 0, "dynamic_VCC"}; %#ok<AGROW>
        end
    end
end

% 7. Grid-Following Converters (GFL, intconv_GFL): CC Dynamic -> Static, PLL Dynamic -> Static
for uCell = {'GFL', 'intconv_GFL'}
    uType = uCell{1};
    ccIdx = find(stateMeta.UnitType == uType & stateMeta.SubGroup == "cc");
    if ~isempty(ccIdx) && all(below(ccIdx))
        curCC = local_get_unit_mode(session, uType, 'dynamic_CC', 'Dynamic');
        if ~strcmpi(curCC, 'Static')
            pendingUnitOverrides.(uType).dynamic_CC = 'Static';
            planRows(end+1, :) = {"converter_cc", uType, string(curCC), "Static", ...
                max(stateRel(ccIdx)), strjoin(stateMeta.Label(ccIdx), ', '), 0, "dynamic_CC"}; %#ok<AGROW>
        end
    end
    pllIdx = find(stateMeta.UnitType == uType & stateMeta.SubGroup == "pll");
    if ~isempty(pllIdx) && all(below(pllIdx))
        curPLL = local_get_unit_mode(session, uType, 'dynamic_PLL', 'Dynamic');
        if ~strcmpi(curPLL, 'Static')
            pendingUnitOverrides.(uType).dynamic_PLL = 'Static';
            planRows(end+1, :) = {"converter_pll", uType, string(curPLL), "Static", ...
                max(stateRel(pllIdx)), strjoin(stateMeta.Label(pllIdx), ', '), 0, "dynamic_PLL"}; %#ok<AGROW>
        end
    end
end

% 8. AC & DC Buses (evaluated last, taking pending line/load/unit choices into account)
[acBusPlans, acBlocked] = local_plan_ac_buses(session, stateMeta, stateRel, below, planRows, pendingUnitOverrides);
planRows    = [planRows; acBusPlans];
blockedMsgs = [blockedMsgs; acBlocked];

[dcBusPlans, dcBlocked] = local_plan_dc_buses(session, stateMeta, stateRel, below, planRows);
planRows    = [planRows; dcBusPlans];
blockedMsgs = [blockedMsgs; dcBlocked];

% Build clean public planTable
if isempty(planRows)
    planTbl = local_empty_plan_table();
else
    planTbl = table( ...
        string(planRows(:, 1)), string(planRows(:, 2)), ...
        string(planRows(:, 3)), string(planRows(:, 4)), ...
        cell2mat(planRows(:, 5)), string(planRows(:, 6)), ...
        'VariableNames', {'Kind','Element','OldComplexity','NewComplexity','MaxRelevance','States'});
end

if dryRun
    out = planTbl;
    return;
end

% Execute planned modifications on session
nApplied = size(planRows, 1);
for i = 1:nApplied
    kindStr = string(planRows{i, 1});
    elemStr = string(planRows{i, 2});
    newComp = string(planRows{i, 4});
    tRow    = planRows{i, 7};
    prmFld  = string(planRows{i, 8});

    switch kindStr
        case {"ac_line", "dc_line"}
            session.tables.Line.Complexity(tRow) = newComp;
            if isfield(session.rawTables, 'Line') && height(session.rawTables.Line) >= tRow && ...
                    ismember('Complexity', session.rawTables.Line.Properties.VariableNames)
                session.rawTables.Line.Complexity = string(session.rawTables.Line.Complexity);
                session.rawTables.Line.Complexity(tRow) = newComp;
            end

        case {"ac_load", "dc_load"}
            session.tables.Bus.Complexity_load(tRow) = newComp;
            if isfield(session.rawTables, 'Bus') && height(session.rawTables.Bus) >= tRow && ...
                    ismember('Complexity_load', session.rawTables.Bus.Properties.VariableNames)
                session.rawTables.Bus.Complexity_load = string(session.rawTables.Bus.Complexity_load);
                session.rawTables.Bus.Complexity_load(tRow) = newComp;
            end

        case {"ac_bus", "dc_bus"}
            session.tables.Bus.Complexity_bus(tRow) = newComp;
            if isfield(session.rawTables, 'Bus') && height(session.rawTables.Bus) >= tRow && ...
                    ismember('Complexity_bus', session.rawTables.Bus.Properties.VariableNames)
                session.rawTables.Bus.Complexity_bus = string(session.rawTables.Bus.Complexity_bus);
                session.rawTables.Bus.Complexity_bus(tRow) = newComp;
            end

        case {"machine", "converter_vcc", "converter_cc", "converter_pll"}
            uName = char(elemStr);
            fName = char(prmFld);
            if ~isfield(session.unitOverrides, uName) || ~isstruct(session.unitOverrides.(uName))
                session.unitOverrides.(uName) = struct();
            end
            session.unitOverrides.(uName).(fName) = char(newComp);
            local_sync_dyr_complexity(session, uName, fName, char(newComp));
    end
end

if nApplied > 0
    session.invalidateDynamicResults();
    session.log('Applied %d complexity reductions from state relevance (threshold = %.3g).', ...
        nApplied, threshold);
else
    session.log('No additional complexity reductions applicable at threshold = %.3g.', threshold);
end

out = struct( ...
    'nSimplified', nApplied, ...
    'threshold', threshold, ...
    'planTable', planTbl, ...
    'blocked', blockedMsgs);
end

% =========================================================================
% Local Helpers
% =========================================================================
function [busPlans, blocked] = local_plan_ac_buses(session, stateMeta, stateRel, below, currentPlans, pendingUnitOv)
busPlans = cell(0, 8);
blocked  = strings(0, 1);

acBusIdx  = find(stateMeta.ElementKind == "ac_bus" & stateMeta.SubGroup == "network");
acBusRows = unique(stateMeta.TableRow(acBusIdx), 'stable');
if isempty(acBusRows), return; end

[TbusAC, TgenAC, TlineAC, ~, ~, ~, Tintconv] = session.splitDomains();

% Determine whether ALL AC lines and AC loads will be static omega after currentPlans
allLinesStatic = true;
if ~isempty(TlineAC) && height(TlineAC) > 0
    for l = 1:height(TlineAC)
        cStr = lower(char(strtrim(TlineAC.Complexity(l))));
        % Check if planned to become static omega
        fb = round(TlineAC.fbus(l)); tb = round(TlineAC.tbus(l));
        eid = sprintf('%d-%d', fb, tb);
        isPlanned = any(strcmp(string(currentPlans(:,1)), "ac_line") & strcmp(string(currentPlans(:,2)), eid));
        if ~isPlanned && ismember(cStr, {'dynamic', 'emt'})
            allLinesStatic = false;
        end
    end
end

allLoadsStatic = true;
for b = 1:height(TbusAC)
    if TbusAC.Pd(b) == 0 && TbusAC.Qd(b) == 0, continue; end
    cStr = lower(char(strtrim(TbusAC.Complexity_load(b))));
    eid  = sprintf('%d', round(TbusAC.bus_i(b)));
    isPlanned = any(strcmp(string(currentPlans(:,1)), "ac_load") & strcmp(string(currentPlans(:,2)), eid));
    if ~isPlanned && contains(cStr, 'dynamic')
        allLoadsStatic = false;
    end
end

% Check candidate AC buses against device algebraic-loop blockers
candRows = [];
for r = acBusRows(:)'
    if r <= 0 || r > height(session.tables.Bus), continue; end
    m = acBusIdx(stateMeta.TableRow(acBusIdx) == r);
    if ~all(below(m)), continue; end
    bId = round(session.tables.Bus.bus_i(r));
    reason = local_ac_device_blocker(session, bId, TgenAC, Tintconv, pendingUnitOv);
    if reason ~= ""
        blocked(end+1, 1) = sprintf('ACBus %d kept dynamic: %s', bId, reason); %#ok<AGROW>
    else
        candRows(end+1) = r; %#ok<AGROW>
    end
end

% In Simulink VFlexP, switching to pure algebraic Ybus ('Static omega' grid
% variant) requires all AC buses, lines, and loads to be static omega; in a
% mixed Dynamic AC Grid, an individual static-omega bus requires local shunt
% admittance and no incident static-omega branch.
allACBusesBecomeStatic = (numel(candRows) == numel(acBusRows)) && allLinesStatic && allLoadsStatic;

for r = candRows
    m    = acBusIdx(stateMeta.TableRow(acBusIdx) == r);
    bId  = round(session.tables.Bus.bus_i(r));
    oldC = char(strtrim(session.tables.Bus.Complexity_bus(r)));
    if ~allACBusesBecomeStatic
        hasShunt = (session.tables.Bus.Bs(r) ~= 0) || (session.tables.Bus.Gs(r) ~= 0);
        hasStatLine = false;
        if ~isempty(TlineAC) && height(TlineAC) > 0
            incL = find(round(TlineAC.fbus) == bId | round(TlineAC.tbus) == bId);
            for il = incL(:)'
                cL = lower(char(strtrim(TlineAC.Complexity(il))));
                eid = sprintf('%d-%d', round(TlineAC.fbus(il)), round(TlineAC.tbus(il)));
                plL = any(strcmp(string(currentPlans(:,1)), "ac_line") & strcmp(string(currentPlans(:,2)), eid));
                if plL || ~ismember(cL, {'dynamic', 'emt'})
                    hasStatLine = true;
                end
            end
        end
        if ~hasShunt || hasStatLine
            blocked(end+1, 1) = sprintf('ACBus %d kept dynamic: mixed EMT/static network requires bus shunt and dynamic incident branches.', bId); %#ok<AGROW>
            continue;
        end
    end
    if ~strcmpi(oldC, 'static omega') && ~strcmpi(oldC, 'static')
        busPlans(end+1, :) = {"ac_bus", sprintf('%d', bId), string(oldC), "static omega", ...
            max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
    end
end
end

function [busPlans, blocked] = local_plan_dc_buses(session, stateMeta, stateRel, below, currentPlans)
busPlans = cell(0, 8);
blocked  = strings(0, 1);

dcBusIdx  = find(stateMeta.ElementKind == "dc_bus" & stateMeta.SubGroup == "network");
dcBusRows = unique(stateMeta.TableRow(dcBusIdx), 'stable');
if isempty(dcBusRows), return; end

[~, ~, ~, TbusDC, ~, TlineDC, Tintconv] = session.splitDomains();

allDCLinesStatic = true;
if ~isempty(TlineDC) && height(TlineDC) > 0
    for l = 1:height(TlineDC)
        cStr = lower(char(strtrim(TlineDC.Complexity(l))));
        eid  = sprintf('%d-%d', round(TlineDC.fbus(l)), round(TlineDC.tbus(l)));
        isPl = any(strcmp(string(currentPlans(:,1)), "dc_line") & strcmp(string(currentPlans(:,2)), eid));
        if ~isPl && contains(cStr, 'dynamic')
            allDCLinesStatic = false;
        end
    end
end

candRows = [];
for r = dcBusRows(:)'
    if r <= 0 || r > height(session.tables.Bus), continue; end
    m = dcBusIdx(stateMeta.TableRow(dcBusIdx) == r);
    if ~all(below(m)), continue; end
    bId = round(session.tables.Bus.bus_i(r));
    hasStation = false;
    if ~isempty(Tintconv) && height(Tintconv) > 0
        hasStation = any(round(Tintconv.fbus) == bId | round(Tintconv.tbus) == bId);
    end
    if hasStation
        blocked(end+1, 1) = sprintf('DCBus %d kept dynamic: interlinking converter DC-link capacitor state required.', bId); %#ok<AGROW>
    else
        candRows(end+1) = r; %#ok<AGROW>
    end
end

allDCBusesStatic = (numel(candRows) == height(TbusDC)) && allDCLinesStatic;
if ~allDCBusesStatic
    for r = candRows
        bId = round(session.tables.Bus.bus_i(r));
        blocked(end+1, 1) = sprintf('DCBus %d kept dynamic: mixed DC network keeps bus capacitance.', bId); %#ok<AGROW>
    end
    return;
end

for r = candRows
    m    = dcBusIdx(stateMeta.TableRow(dcBusIdx) == r);
    bId  = round(session.tables.Bus.bus_i(r));
    oldC = char(strtrim(session.tables.Bus.Complexity_bus(r)));
    if ~strcmpi(oldC, 'DC static')
        busPlans(end+1, :) = {"dc_bus", sprintf('%d', bId), string(oldC), "DC static", ...
            max(stateRel(m)), strjoin(stateMeta.Label(m), ', '), r, ""}; %#ok<AGROW>
    end
end
end

function reason = local_ac_device_blocker(session, busId, TgenAC, Tintconv, pendingUnitOv)
reason = "";
if ~isempty(TgenAC) && height(TgenAC) > 0
    gRows = find(round(TgenAC.bus) == busId);
    for idx = gRows(:)'
        uType = char(TgenAC.type(idx));
        if ismember(uType, {'MGES', 'SG'})
            stMode = local_effective_mode(session, pendingUnitOv, uType, 'dynamic_stator', 'order 8');
            if ~startsWith(lower(stMode), 'order 8')
                reason = "synchronous machine with algebraic stator (order 4) solves current from bus voltage.";
                return;
            end
        elseif strcmp(uType, 'GFL')
            ccMode = local_effective_mode(session, pendingUnitOv, 'GFL', 'dynamic_CC', 'Dynamic');
            if strcmpi(ccMode, 'Static')
                reason = "GFL with Static CC solves current injection from bus voltage.";
                return;
            end
        end
    end
end
if ~isempty(Tintconv) && height(Tintconv) > 0
    cRows = find(round(Tintconv.fbus) == busId | round(Tintconv.tbus) == busId);
    for idx = cRows(:)'
        uType = char(Tintconv.type(idx));
        if strcmp(uType, 'intconv_GFL')
            ccMode = local_effective_mode(session, pendingUnitOv, 'intconv_GFL', 'dynamic_CC', 'Dynamic');
            if strcmpi(ccMode, 'Static')
                reason = "intconv_GFL with Static CC solves current injection from bus voltage.";
                return;
            end
        end
    end
end
end

function modeVal = local_effective_mode(session, pendingOv, uType, fldName, defVal)
if isfield(pendingOv, uType) && isfield(pendingOv.(uType), fldName)
    modeVal = char(string(pendingOv.(uType).(fldName)));
else
    modeVal = local_get_unit_mode(session, uType, fldName, defVal);
end
end

function modeVal = local_get_unit_mode(session, uType, fldName, defVal)
modeVal = defVal;
if isfield(session.unitOverrides, uType) && isfield(session.unitOverrides.(uType), fldName)
    modeVal = char(string(session.unitOverrides.(uType).(fldName)));
    return;
end
if ~isempty(session.dyrData) && isstruct(session.dyrData) && isfield(session.dyrData, 'records')
    recs = session.dyrData.records;
    for k = 1:numel(recs)
        m = recs(k).model;
        c = recs(k).complexity;
        if ismember(uType, {'MGES','SG'}) && ismember(m, {'GENROU','GENSAL'}) && strcmp(fldName, 'dynamic_stator')
            if isfield(c, 'MODE') && strcmpi(c.MODE, 'RMS4'), modeVal = 'order 4'; return; end
            if isfield(c, 'MODE') && strcmpi(c.MODE, 'CLASSICAL2'), modeVal = 'order 2'; return; end
        elseif strcmp(uType, 'GFL') && strcmp(m, 'GFL')
            if strcmp(fldName, 'dynamic_CC') && isfield(c, 'CC'), modeVal = char(c.CC); return; end
            if strcmp(fldName, 'dynamic_PLL') && isfield(c, 'PLL'), modeVal = char(c.PLL); return; end
        elseif strcmp(uType, 'intconv_GFL') && strcmp(m, 'INTCONV_GFL')
            if strcmp(fldName, 'dynamic_CC') && isfield(c, 'CC'), modeVal = char(c.CC); return; end
            if strcmp(fldName, 'dynamic_PLL') && isfield(c, 'PLL'), modeVal = char(c.PLL); return; end
        elseif strcmp(uType, 'VSM') && strcmp(m, 'VSM') && strcmp(fldName, 'dynamic_VCC')
            if isfield(c, 'VCC'), modeVal = char(c.VCC); return; end
        elseif strcmp(uType, 'intconv_VSM') && strcmp(m, 'INTCONV_VSM') && strcmp(fldName, 'dynamic_VCC')
            if isfield(c, 'VCC'), modeVal = char(c.VCC); return; end
        end
    end
end
end

function local_sync_dyr_complexity(session, uType, fldName, newVal)
if isempty(session.dyrData) || ~isstruct(session.dyrData) || ~isfield(session.dyrData, 'records')
    return;
end
for k = 1:numel(session.dyrData.records)
    m = session.dyrData.records(k).model;
    if ismember(uType, {'MGES','SG'}) && ismember(m, {'GENROU','GENSAL'}) && strcmp(fldName, 'dynamic_stator')
        if strcmpi(newVal, 'order 4')
            session.dyrData.records(k).complexity.MODE = 'RMS4';
        elseif strcmpi(newVal, 'order 2')
            session.dyrData.records(k).complexity.MODE = 'CLASSICAL2';
        else
            session.dyrData.records(k).complexity.MODE = 'EMT8';
        end
    elseif ((strcmp(uType, 'GFL') && strcmp(m, 'GFL')) || (strcmp(uType, 'intconv_GFL') && strcmp(m, 'INTCONV_GFL')))
        if strcmp(fldName, 'dynamic_CC')
            session.dyrData.records(k).complexity.CC = upper(newVal);
        elseif strcmp(fldName, 'dynamic_PLL')
            session.dyrData.records(k).complexity.PLL = upper(newVal);
        end
    elseif ((strcmp(uType, 'VSM') && strcmp(m, 'VSM')) || (strcmp(uType, 'intconv_VSM') && strcmp(m, 'INTCONV_VSM')))
        if strcmp(fldName, 'dynamic_VCC')
            session.dyrData.records(k).complexity.VCC = upper(newVal);
        end
    end
end
end

function T = local_empty_plan_table()
T = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', {'Kind','Element','OldComplexity','NewComplexity','MaxRelevance','States'});
end
