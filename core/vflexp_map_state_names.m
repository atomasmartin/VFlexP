function [labels, stateMeta] = vflexp_map_state_names(session, rawStateNames)
%VFLEXP_MAP_STATE_NAMES Map Simulink full state block paths to human-readable labels & metadata.
%
%   [labels, stateMeta] = vflexp_map_state_names(session, rawStateNames)
%
%   Returns:
%     labels    - string(Nx, 1) unique human-readable state names
%     stateMeta - table(Nx, 10) with columns:
%                 Index, Label, Domain, ElementKind, ElementId,
%                 UnitType, TableRow, SubGroup, StateVar, RawPath

rawPaths = string(rawStateNames(:));
nStates  = numel(rawPaths);

if nStates == 0
    labels = strings(0, 1);
    stateMeta = table(zeros(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        strings(0,1), strings(0,1), zeros(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'Index','Label','Domain','ElementKind','ElementId', ...
                          'UnitType','TableRow','SubGroup','StateVar','RawPath'});
    return;
end

[TbusAC, TgenAC, TlineAC, TbusDC, TgenDC, TlineDC, Tintconv] = session.splitDomains();
TBusAll  = session.tables.Bus;
TLineAll = session.tables.Line;
TGenAll  = session.tables.Gen;

% Precompute active network subsets matching vflexp_build_and_init ordering
hasMdlWks = bdIsLoaded('VFlexP');
mdlWks    = [];
if hasMdlWks
    try
        mdlWks = get_param('VFlexP', 'ModelWorkspace');
    catch
        mdlWks = [];
    end
end

busCompAC  = local_ws_or_default(mdlWks, 'bus_complexity', local_ac_bus_codes(TbusAC.Complexity_bus));
lineCompAC = local_ws_or_default(mdlWks, 'line_complexity', local_ac_line_codes(TlineAC.Complexity));
loadIdxAC  = find((TbusAC.Pd ~= 0) | (TbusAC.Qd ~= 0));
if ~isempty(loadIdxAC)
    loadCompAC = local_ws_or_default(mdlWks, 'load_complexity', local_ac_load_codes(TbusAC.Complexity_load(loadIdxAC)));
else
    loadCompAC = [];
end

if ~isempty(TbusDC) && height(TbusDC) > 0
    busCompDC  = local_ws_or_default(mdlWks, 'busDC_complexity', local_dc_bus_codes(TbusDC.Complexity_bus));
    lineCompDC = local_ws_or_default(mdlWks, 'lineDC_complexity', local_dc_line_codes(TlineDC.Complexity));
    loadIdxDC  = find((TbusDC.Pd ~= 0) | (TbusDC.Qd ~= 0));
    if ~isempty(loadIdxDC)
        loadCompDC = local_ws_or_default(mdlWks, 'loadDC_complexity', local_dc_load_codes(TbusDC.Complexity_load(loadIdxDC)));
    else
        loadCompDC = [];
    end
else
    busCompDC  = [];
    lineCompDC = [];
    loadIdxDC  = [];
    loadCompDC = [];
end

% Strip trailing Simulink vector/matrix element suffixes like "(1)" or "(2,1)"
% while preserving explicit 1-based element index when present.
basePaths = strings(nStates, 1);
elemIdx   = nan(nStates, 1);
for k = 1:nStates
    pClean = strrep(char(rawPaths(k)), newline, ' ');
    tokIdx = regexp(pClean, '\((\d+)\)$', 'tokens', 'once');
    if ~isempty(tokIdx)
        elemIdx(k) = str2double(tokIdx{1});
    end
    pBase  = regexprep(pClean, '\(\d+(?:,\s*\d+)*\)$', '');
    basePaths(k) = string(pBase);
end

% Count total occurrences of each base block path (vector/matrix Integrator width)
[uniqPaths, ~, pathGroup] = unique(basePaths, 'stable');
totalPerGroup = accumarray(pathGroup, 1);
seenPerGroup  = zeros(numel(uniqPaths), 1);

labels      = strings(nStates, 1);
domainCol   = strings(nStates, 1);
kindCol     = strings(nStates, 1);
elemIdCol   = strings(nStates, 1);
unitTypeCol = strings(nStates, 1);
tableRowCol = zeros(nStates, 1);
subGrpCol   = strings(nStates, 1);
stateVarCol = strings(nStates, 1);

for k = 1:nStates
    p = char(basePaths(k));
    gIdx = pathGroup(k);
    seenPerGroup(gIdx) = seenPerGroup(gIdx) + 1;
    occ  = seenPerGroup(gIdx);
    nOcc = totalPerGroup(gIdx);

    dom  = "AC";
    kind = "other";
    eid  = "";
    utp  = "";
    trow = 0;
    sgrp = "other";
    svar = local_leaf_name(p);
    lbl  = "";

    % 1. AC Grid Dynamic Buses (dynamic EMT)
    if startsWith(p, 'VFlexP/Grid/AC Grid/Dynamic/Buses/dynamic/Active/Impl/')
        dom  = "AC";
        kind = "ac_bus";
        sgrp = "network";
        dynRows = find(busCompAC == 1);
        if occ <= numel(dynRows)
            bId  = round(TbusAC.bus_i(dynRows(occ)));
            eid  = sprintf('%d', bId);
            trow = local_find_bus_row(TBusAll, bId);
            lbl  = sprintf('ACBus%d.%s', bId, svar);
        else
            lbl  = sprintf('ACBus_dyn(%d).%s', occ, svar);
        end

    % 2. AC Grid Static-Omega Regularization LPF ([3 x nStatOmega] states)
    elseif startsWith(p, 'VFlexP/Grid/AC Grid/Dynamic/Buses/static omega/Active/Impl/')
        dom  = "AC";
        kind = "ac_bus";
        sgrp = "static_reg";
        statRows = find(busCompAC == 3);
        compNames = {'vd_reg', 'vq_reg', 'v0_reg'};
        compIdx = mod(occ - 1, 3) + 1;
        bPos    = floor((occ - 1) / 3) + 1;
        svar    = string(compNames{compIdx});
        if bPos <= numel(statRows)
            bId  = round(TbusAC.bus_i(statRows(bPos)));
            eid  = sprintf('%d', bId);
            trow = local_find_bus_row(TBusAll, bId);
            lbl  = sprintf('ACBus%d.%s', bId, svar);
        else
            lbl  = sprintf('ACBus_reg(%d).%s', bPos, svar);
        end

    % 3. AC Grid Dynamic Lines
    elseif startsWith(p, 'VFlexP/Grid/AC Grid/Dynamic/Lines/dynamic/Active/Impl/')
        dom  = "AC";
        kind = "ac_line";
        sgrp = "network";
        dynL = find(lineCompAC == 1);
        if occ <= numel(dynL)
            lr   = dynL(occ);
            fb   = round(TlineAC.fbus(lr));
            tb   = round(TlineAC.tbus(lr));
            eid  = sprintf('%d-%d', fb, tb);
            trow = local_find_line_row(TLineAll, fb, tb, false);
            lbl  = sprintf('ACLine(%s).%s', eid, svar);
        else
            lbl  = sprintf('ACLine_dyn(%d).%s', occ, svar);
        end

    % 4. AC Grid Dynamic Loads
    elseif startsWith(p, 'VFlexP/Grid/AC Grid/Dynamic/Loads/')
        dom  = "AC";
        kind = "ac_load";
        sgrp = "network";
        relLd = extractAfter(p, 'VFlexP/Grid/AC Grid/Dynamic/Loads/');
        slashPos = find(relLd == '/', 1);
        varName  = lower(relLd(1:slashPos-1));
        ldCode = 0;
        if strcmp(varName, 'series rl dynamic'),   ldCode = 1;
        elseif strcmp(varName, 'dynamic pq'),      ldCode = 5;
        elseif strcmp(varName, 'parallel rl dynamic'), ldCode = 6;
        end
        ldSub = loadIdxAC(loadCompAC == ldCode);
        if ~isempty(ldSub) && nOcc >= numel(ldSub)
            mult = max(1, round(nOcc / numel(ldSub)));
            uPos = floor((occ - 1) / mult) + 1;
            cPos = mod(occ - 1, mult) + 1;
            if mult > 1
                svar = sprintf('%s_%d', svar, cPos);
            end
            if uPos <= numel(ldSub)
                bId  = round(TbusAC.bus_i(ldSub(uPos)));
                eid  = sprintf('%d', bId);
                trow = local_find_bus_row(TBusAll, bId);
                lbl  = sprintf('ACLoad@Bus%d.%s', bId, svar);
            end
        end
        if lbl == ""
            lbl = sprintf('ACLoad(%d).%s', occ, svar);
        end

    % 5. DC Grid Dynamic Buses
    elseif startsWith(p, 'VFlexP/Grid/DC Grid/Dynamic/Buses/dynamic/Active/Impl/')
        dom  = "DC";
        kind = "dc_bus";
        sgrp = "network";
        svar = "v_dc";
        dynB = find(busCompDC == 1);
        if occ <= numel(dynB)
            bId  = round(TbusDC.bus_i(dynB(occ)));
            eid  = sprintf('%d', bId);
            trow = local_find_bus_row(TBusAll, bId);
            lbl  = sprintf('DCBus%d.v_dc', bId);
        else
            lbl  = sprintf('DCBus_dyn(%d).v_dc', occ);
        end

    % 6. DC Grid Dynamic Lines
    elseif startsWith(p, 'VFlexP/Grid/DC Grid/Dynamic/Lines/Dynamic/Active/Impl/')
        dom  = "DC";
        kind = "dc_line";
        sgrp = "network";
        svar = "i_dc";
        dynL = find(lineCompDC == 1);
        if occ <= numel(dynL)
            lr   = dynL(occ);
            fb   = round(TlineDC.fbus(lr));
            tb   = round(TlineDC.tbus(lr));
            eid  = sprintf('%d-%d', fb, tb);
            trow = local_find_line_row(TLineAll, fb, tb, true);
            lbl  = sprintf('DCLine(%s).i_dc', eid);
        else
            lbl  = sprintf('DCLine_dyn(%d).i_dc', occ);
        end

    % 7. DC Grid DynamicP Loads
    elseif startsWith(p, 'VFlexP/Grid/DC Grid/Dynamic/Loads/DynamicP/Active/Impl/')
        dom  = "DC";
        kind = "dc_load";
        sgrp = "network";
        svar = "PI_load";
        ldSub = loadIdxDC(loadCompDC == 4);
        if occ <= numel(ldSub)
            bId  = round(TbusDC.bus_i(ldSub(occ)));
            eid  = sprintf('%d', bId);
            trow = local_find_bus_row(TBusAll, bId);
            lbl  = sprintf('DCLoad@Bus%d.PI', bId);
        else
            lbl  = sprintf('DCLoad(%d).PI', occ);
        end

    % 8. AC Generators (Slot_<U>/Active/<U>/<subPath>)
    elseif startsWith(p, 'VFlexP/Grid/AC Generators/Slot_')
        dom  = "AC";
        kind = "ac_gen";
        [utp, subPath] = local_parse_slot_path(p, 'VFlexP/Grid/AC Generators/Slot_');
        rowsU = find(strcmp(TgenAC.type, utp));
        nU    = max(1, numel(rowsU));
        if strcmp(utp, "MGES") && endsWith(subPath, 'AVR_LeadLag/x_f') && ~isempty(mdlWks)
            subU = local_mges_leadlag_units(mdlWks, nU);
            mult = 1;
            if ~isnan(elemIdx(k)) && elemIdx(k) >= 1 && elemIdx(k) <= numel(subU)
                uPos = subU(elemIdx(k));
            elseif occ <= numel(subU)
                uPos = subU(occ);
            else
                uPos = min(nU, occ);
            end
            cPos = 1;
        elseif nOcc < nU && ~isnan(elemIdx(k)) && elemIdx(k) >= 1 && elemIdx(k) <= nU
            mult = 1;
            uPos = elemIdx(k);
            cPos = 1;
        else
            mult  = max(1, round(nOcc / nU));
            uPos  = min(nU, floor((occ - 1) / mult) + 1);
            cPos  = mod(occ - 1, mult) + 1;
        end
        [sgrp, svar] = local_classify_unit_subpath(utp, subPath, mult, cPos);
        if ~isempty(rowsU)
            bId  = round(TgenAC.bus(rowsU(uPos)));
            eid  = sprintf('%d', bId);
            trow = local_find_gen_row(TGenAll, bId, utp, uPos, rowsU, TgenAC);
            lbl  = sprintf('%s@Bus%d.%s', utp, bId, svar);
        else
            lbl  = sprintf('%s(%d).%s', utp, uPos, svar);
        end

    % 9. DC Generators (Slot_<U>/Active/<U>/<subPath>)
    elseif startsWith(p, 'VFlexP/Grid/DC Generators/Slot_')
        dom  = "DC";
        kind = "dc_gen";
        [utp, subPath] = local_parse_slot_path(p, 'VFlexP/Grid/DC Generators/Slot_');
        rowsD = find(strcmp(TgenDC.type, utp));
        nD    = max(1, numel(rowsD));
        mult  = max(1, round(nOcc / nD));
        uPos  = min(nD, floor((occ - 1) / mult) + 1);
        cPos  = mod(occ - 1, mult) + 1;
        [sgrp, svar] = local_classify_unit_subpath(utp, subPath, mult, cPos);
        if ~isempty(rowsD)
            bId  = round(TgenDC.bus(rowsD(uPos)));
            eid  = sprintf('%d', bId);
            trow = local_find_gen_row(TGenAll, bId, utp, uPos, rowsD, TgenDC);
            lbl  = sprintf('%s@Bus%d.%s', utp, bId, svar);
        else
            lbl  = sprintf('%s(%d).%s', utp, uPos, svar);
        end

    % 10. DC/AC Interlinking Converters (Slot_<U>/Active/<U>/<subPath>)
    elseif startsWith(p, 'VFlexP/Grid/DCAC Converters/Slot_')
        dom  = "Hybrid";
        kind = "intconv";
        [utp, subPath] = local_parse_slot_path(p, 'VFlexP/Grid/DCAC Converters/Slot_');
        rowsIC = find(strcmp(Tintconv.type, utp));
        nC     = max(1, numel(rowsIC));
        if strcmp(utp, "intconv_GFL") && endsWith(subPath, 'VdcCtrl/Active/I')
            subC = find(round(Tintconv.CONV(rowsIC)) == 2);
            mult = 1;
            if ~isnan(elemIdx(k)) && elemIdx(k) >= 1 && elemIdx(k) <= numel(subC)
                uPos = subC(elemIdx(k));
            elseif occ <= numel(subC)
                uPos = subC(occ);
            else
                uPos = min(nC, occ);
            end
            cPos = 1;
        else
            mult   = max(1, round(nOcc / nC));
            uPos   = min(nC, floor((occ - 1) / mult) + 1);
            cPos   = mod(occ - 1, mult) + 1;
        end
        [sgrp, svar] = local_classify_unit_subpath(utp, subPath, mult, cPos);
        if ~isempty(rowsIC)
            icRow = rowsIC(uPos);
            fb    = round(Tintconv.fbus(icRow));
            tb    = round(Tintconv.tbus(icRow));
            eid   = sprintf('%d-%d', fb, tb);
            trow  = icRow;
            lbl   = sprintf('%s(%s).%s', utp, eid, svar);
        else
            lbl   = sprintf('%s(%d).%s', utp, uPos, svar);
        end

    % 11. Secondary / Tertiary Controls & Fallback
    elseif contains(p, '/Secondary Control/')
        dom  = "Control";
        kind = "sec_ctrl";
        sgrp = "sec_ctrl";
        lbl  = sprintf('SecCtrl.%s(%d)', svar, occ);
    elseif contains(p, '/Tertiary Control/')
        dom  = "Control";
        kind = "ter_ctrl";
        sgrp = "ter_ctrl";
        lbl  = sprintf('TerCtrl.%s(%d)', svar, occ);
    else
        shortP = erase(p, 'VFlexP/Grid/');
        lbl = sprintf('%s(%d)', strrep(shortP, '/', '.'), occ);
    end

    labels(k)      = lbl;
    domainCol(k)   = dom;
    kindCol(k)     = kind;
    elemIdCol(k)   = eid;
    unitTypeCol(k) = utp;
    tableRowCol(k) = trow;
    subGrpCol(k)   = sgrp;
    stateVarCol(k) = svar;
end

% Guarantee strict uniqueness of labels when multiple units share a bus & type
labels = local_ensure_unique_labels(labels);

stateMeta = table((1:nStates)', labels, domainCol, kindCol, elemIdCol, ...
    unitTypeCol, tableRowCol, subGrpCol, stateVarCol, rawPaths, ...
    'VariableNames', {'Index','Label','Domain','ElementKind','ElementId', ...
                      'UnitType','TableRow','SubGroup','StateVar','RawPath'});
end

% =========================================================================
% Local Helpers
% =========================================================================
function v = local_ws_or_default(mdlWks, varName, defVal)
v = defVal;
if isempty(mdlWks), return; end
try
    cand = mdlWks.getVariable(varName);
    if isnumeric(cand) && numel(cand) == numel(defVal)
        v = reshape(double(cand), 1, []);
    end
catch
end
end

function subU = local_mges_leadlag_units(mdlWks, nU)
subU = 1:nU;
try
    prm = mdlWks.getVariable('prm_MGES');
    if isstruct(prm) && isfield(prm, 'T_1')
        hit = find(double(prm.T_1) ~= 0);
        if ~isempty(hit)
            subU = reshape(hit, 1, []);
        end
    end
catch
end
end

function leaf = local_leaf_name(p)
parts = strsplit(p, '/');
tok   = erase(strtrim(parts{end}), '\');
leaf  = string(regexprep(tok, '[^A-Za-z0-9_]', '_'));
end

function [uType, subPath] = local_parse_slot_path(p, prefix)
rest   = extractAfter(p, prefix);
parts  = strsplit(rest, '/');
uType  = string(parts{1});
% Pattern: <U>/Active/<U>/<subPath...>
if numel(parts) >= 4 && strcmp(parts{2}, 'Active')
    subPath = strjoin(parts(4:end), '/');
else
    subPath = strjoin(parts(2:end), '/');
end
end

function [subGrp, stateVar] = local_classify_unit_subpath(~, subPath, mult, cPos)
sp = char(subPath);
parts = strsplit(sp, '/');
leaf  = regexprep(erase(parts{end}, '\'), '[^A-Za-z0-9_]', '_');

% Clean common wrapper tokens from subPath to form a readable stateVar
cleanParts = parts(~ismember(parts, {'Active', 'Impl', 'Stator circuit', 'dynamic', 'Rotor axis'}));
if numel(cleanParts) >= 2
    prevTok  = regexprep(erase(cleanParts{end-1}, '\'), '[^A-Za-z0-9_]', '_');
    stateVar = string(sprintf('%s_%s', prevTok, leaf));
else
    stateVar = string(leaf);
end
if mult > 1
    stateVar = sprintf('%s_%d', stateVar, cPos);
end

subGrp = "control";
if startsWith(sp, 'Electromagnetic/')
    if contains(sp, 'order 8') && ismember(leaf, {'psis_d','psis_q','psikd2','psikq3','Integrator'})
        subGrp = "stator_emt";
        stateVar = string(leaf);
    elseif ismember(leaf, {'psifd','psikd1','psikq1','psikq2'})
        subGrp = "rotor_flux";
        stateVar = string(leaf);
    elseif contains(sp, 'order 2')
        subGrp = "stator_emt";
        stateVar = string(leaf);
    else
        subGrp = "electromagnetic";
    end
elseif startsWith(sp, 'Rotor axis/') || contains(sp, 'omega_r') || strcmp(leaf, 'theta') || strcmp(leaf, 'frame')
    subGrp = "swing";
elseif startsWith(sp, 'VCC/')
    subGrp = "vcc";
elseif startsWith(sp, 'CC/')
    subGrp = "cc";
elseif startsWith(sp, 'PLL/')
    subGrp = "pll";
elseif startsWith(sp, 'RL_out/') || startsWith(sp, 'RL out/')
    subGrp = "rl_out";
elseif startsWith(sp, 'AVR') || contains(sp, 'LeadLag') || contains(sp, 'Exciter')
    subGrp = "avr";
elseif startsWith(sp, 'Gov') || contains(sp, 'Turbine') || contains(sp, 'Servo')
    subGrp = "gov";
elseif startsWith(sp, 'PSS') || contains(sp, 'Washout')
    subGrp = "pss";
elseif startsWith(sp, 'PC/')
    subGrp = "power_ctrl";
elseif contains(sp, 'LPF') || contains(sp, 'Filter')
    subGrp = "filter";
end
end

function labels = local_ensure_unique_labels(labels)
[u, ~, ic] = unique(labels, 'stable');
counts = accumarray(ic, 1);
seen   = zeros(numel(u), 1);
for k = 1:numel(labels)
    g = ic(k);
    if counts(g) > 1
        seen(g) = seen(g) + 1;
        labels(k) = sprintf('%s#%d', labels(k), seen(g));
    end
end
end

function r = local_find_bus_row(TBusAll, busId)
r = find(round(TBusAll.bus_i) == busId, 1);
if isempty(r), r = 0; end
end

function r = local_find_line_row(TLineAll, fb, tb, isDC)
if isempty(TLineAll) || height(TLineAll) == 0
    r = 0; return;
end
dcMask = startsWith(strtrim(TLineAll.Complexity), "DC", 'IgnoreCase', true);
if isDC
    cand = find(dcMask & round(TLineAll.fbus) == fb & round(TLineAll.tbus) == tb, 1);
else
    cand = find(~dcMask & round(TLineAll.fbus) == fb & round(TLineAll.tbus) == tb, 1);
end
if isempty(cand), r = 0; else, r = cand; end
end

function r = local_find_gen_row(TGenAll, busId, uType, uPos, rowsSub, TgenSub)
% Match exact row in session.tables.Gen
subRow = rowsSub(uPos);
cand   = find(round(TGenAll.bus) == busId & strcmp(TGenAll.type, uType));
if isempty(cand)
    r = 0;
elseif isscalar(cand)
    r = cand;
else
    sameBefore = sum(round(TgenSub.bus(1:subRow)) == busId & strcmp(TgenSub.type(1:subRow), uType));
    r = cand(min(sameBefore, numel(cand)));
end
end

function c = local_ac_bus_codes(strs)
n = numel(strs); c = 3 * ones(1, n);
opts = {'dynamic', 'static', 'static omega', 'only r'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if strcmp(s, 'emt'), s = 'dynamic'; elseif strcmp(s, 'rms'), s = 'static omega'; end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_ac_line_codes(strs)
n = numel(strs); c = 3 * ones(1, n);
opts = {'dynamic', 'static', 'static omega'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if strcmp(s, 'emt'), s = 'dynamic'; elseif strcmp(s, 'rms'), s = 'static omega'; end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_ac_load_codes(strs)
n = numel(strs); c = 3 * ones(1, n);
opts = {'series rl dynamic', 'series rl static', 'series rl static omega', ...
        'pq', 'dynamic pq', 'parallel rl dynamic', 'parallel rl static', ...
        'parallel rl static omega', 'constant current'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'emt','dynamic'}), s = 'series rl dynamic';
    elseif ismember(s, {'rms','static omega'}), s = 'series rl static omega';
    elseif strcmp(s, 'static'), s = 'series rl static';
    end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_dc_bus_codes(strs)
n = numel(strs); c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc dynamic', 'dc emt', 'dynamic', 'emt'}), c(k) = 1; end
end
end

function c = local_dc_line_codes(strs)
n = numel(strs); c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc dynamic', 'dc emt', 'dynamic', 'emt'}), c(k) = 1; end
end
end

function c = local_dc_load_codes(strs)
n = numel(strs); c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc constr', 'constr'}), c(k) = 1;
    elseif ismember(s, {'dc constp', 'constp', 'dc static', 'dc static omega', 'dc rms'}), c(k) = 2;
    elseif ismember(s, {'dc consti', 'consti'}), c(k) = 3;
    elseif ismember(s, {'dc dynamicp', 'dynamicp', 'dc dynamic', 'dc emt'}), c(k) = 4;
    end
end
end
