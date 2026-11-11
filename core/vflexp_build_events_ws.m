function evWs = vflexp_build_events_ws(session)
%VFLEXP_BUILD_EVENTS_WS Build FromWorkspace structs and t_close_lines from session.events.
%
%   evWs = vflexp_build_events_ws(session)
%   Returns struct with fields:
%     ev_bus_dP, ev_bus_dQ, ev_dc_dI, ev_gen_dref, t_close_lines

[TbusAC, TgenAC, TlineAC, TbusDC, TgenDC, ~, Tintconv, ~, Sb] = session.splitDomains();

nAC  = max(1, height(TbusAC));
nDC  = max(1, height(TbusDC));
nAg  = max(1, height(TgenAC) + height(TgenDC) + height(Tintconv));
nLin = max(1, height(TlineAC));

tEnd = 10.0;
if isfield(session.config, 'sim') && isfield(session.config.sim, 't_end') && ~isempty(session.config.sim.t_end)
    tEnd = max(1e-3, double(session.config.sim.t_end));
end

t_close_lines = zeros(1, nLin);
evTbl = session.eventList;

if isempty(evTbl) || height(evTbl) == 0
    evWs = local_zero_structs(nAC, nDC, nAg, t_close_lines, tEnd);
    return;
end

% Filter enabled rows
enMask = logical(evTbl.Enabled);
evTbl = evTbl(enMask, :);
if height(evTbl) == 0
    evWs = local_zero_structs(nAC, nDC, nAg, t_close_lines, tEnd);
    return;
end

% Handle AC Line closing events first
isLineEv = strcmpi(strtrim(string(evTbl.TargetType)), "AC Line") | ...
           strcmpi(strtrim(string(evTbl.Channel)), "close");
for r = find(isLineEv(:))'
    lIdx = local_resolve_line(TlineAC, evTbl.TargetId(r));
    if ~isempty(lIdx)
        t_close_lines(lIdx) = max(0, double(evTbl.Time(r)));
    end
end

sigTbl = evTbl(~isLineEv, :);
if height(sigTbl) == 0
    evWs = local_zero_structs(nAC, nDC, nAg, t_close_lines, tEnd);
    return;
end

% Collect switching instants with sharp +/- eps edges for piecewise-constant steps/pulses
epsT = 1e-7;
pts = [0.0; tEnd];
for r = 1:height(sigTbl)
    t0 = max(0.0, double(sigTbl.Time(r)));
    dur = double(sigTbl.Duration(r));
    if t0 > 2*epsT
        pts = [pts; t0 - epsT; t0]; %#ok<AGROW>
    else
        pts = [pts; t0]; %#ok<AGROW>
    end
    if isfinite(dur) && dur > 0
        t1 = t0 + dur;
        pts = [pts; max(0, t1 - epsT); t1]; %#ok<AGROW>
    end
end
tVec = unique(round(pts, 9));
tVec = sort(tVec(:));
% Ensure strictly increasing time vector
dt = diff(tVec);
tVec([false; dt <= 0]) = [];
nT = numel(tVec);

valBusP = zeros(1, nAC, nT);
valBusQ = zeros(1, nAC, nT);
valDcI  = zeros(1, nDC, nT);
valRef  = zeros(6, nAg, nT);
ev_bus_constZ = 0;

nGenAC = height(TgenAC);
nGenDC = height(TgenDC);
nIC    = height(Tintconv);

isVfRow = strcmpi(strtrim(string(sigTbl.TargetType)), "vflexpy");
legTbl  = sigTbl(~isVfRow, :);
vfTbl   = sigTbl(isVfRow, :);

for r = 1:height(legTbl)
    t0  = max(0.0, double(legTbl.Time(r)));
    dur = double(legTbl.Duration(r));
    mag = double(legTbl.Magnitude(r));
    tt  = lower(strtrim(char(string(legTbl.TargetType(r)))));
    ch  = lower(strtrim(char(string(legTbl.Channel(r)))));
    tid = strtrim(char(string(legTbl.TargetId(r))));

    if isfinite(dur) && dur > 0
        activeT = (tVec >= t0) & (tVec < t0 + dur - 0.5*epsT);
    else
        activeT = (tVec >= t0);
    end
    if ~any(activeT)
        continue;
    end

    switch tt
        case {'ac bus', 'bus'}
            bIdx = local_resolve_bus(TbusAC, tid);
            if ~isempty(bIdx)
                if contains(ch, '_z') || contains(ch, 'constz')
                    ev_bus_constZ = 1;
                end
                if contains(ch, 'dp_load') || strcmp(ch, 'pload_mw')
                    valBusP(1, bIdx, activeT) = valBusP(1, bIdx, activeT) - (mag / Sb);
                elseif contains(ch, 'dq_load') || strcmp(ch, 'qload_mvar')
                    valBusQ(1, bIdx, activeT) = valBusQ(1, bIdx, activeT) - (mag / Sb);
                elseif contains(ch, 'dq')
                    valBusQ(1, bIdx, activeT) = valBusQ(1, bIdx, activeT) + mag;
                else
                    valBusP(1, bIdx, activeT) = valBusP(1, bIdx, activeT) + mag;
                end
            end

        case 'dc bus'
            bIdx = local_resolve_bus(TbusDC, tid);
            if ~isempty(bIdx)
                vdc0 = max(TbusDC.Vm(bIdx), 0.1);
                if contains(ch, 'dp_load') || strcmp(ch, 'pload_mw')
                    valDcI(1, bIdx, activeT) = valDcI(1, bIdx, activeT) - (mag / (Sb * vdc0));
                else
                    valDcI(1, bIdx, activeT) = valDcI(1, bIdx, activeT) + mag;
                end
            end

        case {'ac generator', 'generator', 'gen'}
            gIdx = local_resolve_gen(TgenAC, tid);
            if ~isempty(gIdx)
                rCh = local_ref_row(ch);
                valRef(rCh, gIdx, activeT) = valRef(rCh, gIdx, activeT) + mag;
            end

        case 'dc generator'
            gIdx = local_resolve_gen(TgenDC, tid);
            if ~isempty(gIdx)
                agIdx = nGenAC + gIdx;
                rCh = local_ref_row(ch);
                valRef(rCh, agIdx, activeT) = valRef(rCh, agIdx, activeT) + mag;
            end

        case {'converter', 'intconv'}
            icIdx = local_resolve_line(Tintconv, tid);
            if ~isempty(icIdx)
                agIdx = nGenAC + nGenDC + icIdx;
                rCh = local_ref_row(ch);
                valRef(rCh, agIdx, activeT) = valRef(rCh, agIdx, activeT) + mag;
            end
    end
end

if height(vfTbl) > 0
    [uLoadP0_ac, uLoadQ0_ac, uLoadP0_dc, uAgRef0] = local_steady_inputs( ...
        session, TbusAC, TgenAC, TbusDC, TgenDC, Tintconv, Sb, nAC, nDC, nAg);
    uLoadP_ac = repmat(uLoadP0_ac, 1, nT);
    uLoadQ_ac = repmat(uLoadQ0_ac, 1, nT);
    uLoadP_dc = repmat(uLoadP0_dc, 1, nT);
    uAgRef    = repmat(uAgRef0, 1, 1, nT);

    [~, ord] = sortrows([double(vfTbl.Time(:)), (1:height(vfTbl))']);
    for idxOrd = 1:numel(ord)
        r   = ord(idxOrd);
        t0  = max(0.0, double(vfTbl.Time(r)));
        dur = double(vfTbl.Duration(r));
        mag = double(vfTbl.Magnitude(r));
        op  = strtrim(char(string(vfTbl.Channel(r))));
        tid = strtrim(char(string(vfTbl.TargetId(r))));

        if isfinite(dur) && dur > 0
            activeT = (tVec >= t0) & (tVec < t0 + dur - 0.5*epsT);
        else
            activeT = (tVec >= t0);
        end
        if ~any(activeT)
            continue;
        end

        tok = regexp(tid, '^([A-Za-z_]+)\.([^.]+)\.([A-Za-z0-9_]+)$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        kind    = lower(tok{1});
        sel     = strtrim(tok{2});
        varName = lower(tok{3});

        switch kind
            case 'load'
                if ~isempty(TbusAC) && height(TbusAC) > 0
                    if strcmp(sel, '*')
                        bIdxAC = find((TbusAC.Pd ~= 0) | (TbusAC.Qd ~= 0));
                    else
                        bIdxAC = local_resolve_bus_exact(TbusAC, sel);
                    end
                else
                    bIdxAC = [];
                end
                if ~isempty(bIdxAC)
                    ev_bus_constZ = 1;
                    if strcmp(varName, 'p')
                        uLoadP_ac(bIdxAC, activeT) = local_apply_op(uLoadP_ac(bIdxAC, activeT), op, mag);
                    elseif strcmp(varName, 'q')
                        uLoadQ_ac(bIdxAC, activeT) = local_apply_op(uLoadQ_ac(bIdxAC, activeT), op, mag);
                    end
                end

                if ~isempty(TbusDC) && height(TbusDC) > 0 && strcmp(varName, 'p')
                    if strcmp(sel, '*')
                        bIdxDC = find(TbusDC.Pd ~= 0);
                    else
                        bIdxDC = local_resolve_bus_exact(TbusDC, sel);
                    end
                    if ~isempty(bIdxDC)
                        uLoadP_dc(bIdxDC, activeT) = local_apply_op(uLoadP_dc(bIdxDC, activeT), op, mag);
                    end
                end

            case 'gen'
                rCh = local_vf_ref_row(varName);
                agTargets = [];
                if nGenAC > 0
                    if strcmp(sel, '*')
                        agTargets = [agTargets, 1:nGenAC]; %#ok<AGROW>
                    else
                        gIdx = local_resolve_gen_exact(TgenAC, sel);
                        agTargets = [agTargets, reshape(gIdx, 1, [])]; %#ok<AGROW>
                    end
                end
                if nGenDC > 0
                    if strcmp(sel, '*')
                        agTargets = [agTargets, nGenAC + (1:nGenDC)]; %#ok<AGROW>
                    else
                        gIdx = local_resolve_gen_exact(TgenDC, sel);
                        agTargets = [agTargets, nGenAC + reshape(gIdx, 1, [])]; %#ok<AGROW>
                    end
                end
                if ~isempty(agTargets)
                    uAgRef(rCh, agTargets, activeT) = local_apply_op(uAgRef(rCh, agTargets, activeT), op, mag);
                end

            case {'conv', 'intconv'}
                if nIC > 0
                    rCh = local_vf_ref_row(varName);
                    if strcmp(sel, '*')
                        icIdx = 1:nIC;
                    else
                        icIdx = local_resolve_line(Tintconv, sel);
                    end
                    if ~isempty(icIdx)
                        agTargets = nGenAC + nGenDC + reshape(icIdx, 1, []);
                        uAgRef(rCh, agTargets, activeT) = local_apply_op(uAgRef(rCh, agTargets, activeT), op, mag);
                    end
                end
        end
    end

    dLoadP_ac = uLoadP_ac - repmat(uLoadP0_ac, 1, nT);
    dLoadQ_ac = uLoadQ_ac - repmat(uLoadQ0_ac, 1, nT);
    valBusP   = valBusP - reshape(dLoadP_ac, 1, nAC, nT);
    valBusQ   = valBusQ - reshape(dLoadQ_ac, 1, nAC, nT);

    if ~isempty(TbusDC) && height(TbusDC) > 0
        vdc0 = max(reshape(double(TbusDC.Vm), nDC, 1), 0.1);
        dLoadI_dc = (uLoadP_dc - repmat(uLoadP0_dc, 1, nT)) ./ repmat(vdc0, 1, nT);
        valDcI    = valDcI - reshape(dLoadI_dc, 1, nDC, nT);
    end

    valRef = valRef + (uAgRef - repmat(uAgRef0, 1, 1, nT));
end

evWs = struct();
evWs.ev_bus_dP   = struct('time', tVec, 'signals', struct('values', valBusP, 'dimensions', [1, nAC]));
evWs.ev_bus_dQ   = struct('time', tVec, 'signals', struct('values', valBusQ, 'dimensions', [1, nAC]));
evWs.ev_dc_dI    = struct('time', tVec, 'signals', struct('values', valDcI,  'dimensions', [1, nDC]));
evWs.ev_gen_dref = struct('time', tVec, 'signals', struct('values', valRef,  'dimensions', [6, nAg]));
evWs.ev_bus_constZ = ev_bus_constZ;
evWs.t_close_lines = t_close_lines;
end

% -------------------------------------------------------------------------
function evWs = local_zero_structs(nAC, nDC, nAg, t_close_lines, tEnd)
tVec = [0.0; max(1e-3, tEnd)];
evWs = struct();
evWs.ev_bus_dP   = struct('time', tVec, 'signals', struct('values', zeros(1, nAC, 2), 'dimensions', [1, nAC]));
evWs.ev_bus_dQ   = struct('time', tVec, 'signals', struct('values', zeros(1, nAC, 2), 'dimensions', [1, nAC]));
evWs.ev_dc_dI    = struct('time', tVec, 'signals', struct('values', zeros(1, nDC, 2), 'dimensions', [1, nDC]));
evWs.ev_gen_dref = struct('time', tVec, 'signals', struct('values', zeros(6, nAg, 2), 'dimensions', [6, nAg]));
evWs.ev_bus_constZ = 0;
evWs.t_close_lines = t_close_lines;
end

function out = local_apply_op(cur, op, val)
switch op
    case '='
        out = repmat(val, size(cur));
    case '+='
        out = cur + val;
    case '-='
        out = cur - val;
    case '*='
        out = cur * val;
    otherwise
        out = cur;
end
end

function [uLoadP0_ac, uLoadQ0_ac, uLoadP0_dc, uAgRef0] = local_steady_inputs( ...
    session, TbusAC, TgenAC, TbusDC, TgenDC, Tintconv, Sb, nAC, nDC, nAg)
uLoadP0_ac = zeros(nAC, 1);
uLoadQ0_ac = zeros(nAC, 1);
if ~isempty(TbusAC) && height(TbusAC) > 0
    uLoadP0_ac = double(TbusAC.Pd(:)) / Sb;
    uLoadQ0_ac = double(TbusAC.Qd(:)) / Sb;
end

uLoadP0_dc = zeros(nDC, 1);
if ~isempty(TbusDC) && height(TbusDC) > 0
    uLoadP0_dc = double(TbusDC.Pd(:)) / Sb;
end

uAgRef0 = ones(6, nAg);
uAgRef0(2, :) = 0;
uAgRef0(4, :) = 0;

nGenAC = height(TgenAC);
nGenDC = height(TgenDC);
nIC    = height(Tintconv);

for g = 1:nGenAC
    uAgRef0(2, g) = double(TgenAC.Pg(g)) / Sb;
    uAgRef0(4, g) = double(TgenAC.Qg(g)) / Sb;
    bPos = find(TbusAC.bus_i == TgenAC.bus(g), 1);
    if ~isempty(bPos)
        uAgRef0(3, g) = double(TbusAC.Vm(bPos));
    end
end
for g = 1:nGenDC
    ag = nGenAC + g;
    uAgRef0(2, ag) = double(TgenDC.Pg(g)) / Sb;
    if ~isempty(TbusDC) && height(TbusDC) > 0
        bPos = find(TbusDC.bus_i == TgenDC.bus(g), 1);
        if ~isempty(bPos)
            uAgRef0(3, ag) = double(TbusDC.Vm(bPos));
        end
    end
end
if nIC > 0
    hasICPF = ~isempty(session.pfResult) && isstruct(session.pfResult) && ...
              isfield(session.pfResult, 'intconv_PF') && ...
              isfield(session.pfResult.intconv_PF, 'include') && ...
              logical(session.pfResult.intconv_PF.include);
    for c = 1:nIC
        ag = nGenAC + nGenDC + c;
        if hasICPF
            uAgRef0(2, ag) = -double(session.pfResult.intconv_PF.P(c));
            uAgRef0(3, ag) =  double(session.pfResult.intconv_PF.vmod(c));
            uAgRef0(4, ag) = -double(session.pfResult.intconv_PF.Q(c));
        else
            uAgRef0(2, ag) = -double(Tintconv.PT(c)) / Sb;
            uAgRef0(4, ag) = -double(Tintconv.QT(c)) / Sb;
            bPos = find(TbusAC.bus_i == Tintconv.tbus(c), 1);
            if ~isempty(bPos)
                uAgRef0(3, ag) = double(TbusAC.Vm(bPos));
            end
        end
    end
end
end

function bIdx = local_resolve_bus_exact(TBus, tid)
bIdx = [];
if isempty(TBus) || height(TBus) == 0, return; end
num = str2double(tid);
if ~isnan(num)
    bIdx = find(TBus.bus_i == num);
end
end

function gIdx = local_resolve_gen_exact(TGen, tid)
gIdx = [];
if isempty(TGen) || height(TGen) == 0, return; end
num = str2double(tid);
if ~isnan(num)
    gIdx = find(TGen.bus == num);
end
end

function rCh = local_vf_ref_row(varName)
switch varName
    case {'omega_ref', 'w_ref', 'omega'}
        rCh = 1;
    case {'v_ref', 'v'}
        rCh = 3;
    case {'q_ref', 'q'}
        rCh = 4;
    otherwise
        rCh = 2;
end
end

function bIdx = local_resolve_bus(TBus, tid)
bIdx = [];
if isempty(TBus) || height(TBus) == 0, return; end
num = str2double(tid);
if ~isnan(num)
    bIdx = find(TBus.bus_i == num, 1);
    if isempty(bIdx) && num >= 1 && num <= height(TBus)
        bIdx = round(num);
    end
end
end

function gIdx = local_resolve_gen(TGen, tid)
gIdx = [];
if isempty(TGen) || height(TGen) == 0, return; end
% Support "bus:id", bus number, or "#row"
if startsWith(tid, '#')
    r = str2double(tid(2:end));
    if ~isnan(r) && r >= 1 && r <= height(TGen)
        gIdx = round(r);
    end
    return;
end
num = str2double(tid);
if ~isnan(num)
    gIdx = find(TGen.bus == num, 1);
    if isempty(gIdx) && num >= 1 && num <= height(TGen)
        gIdx = round(num);
    end
end
end

function lIdx = local_resolve_line(TLine, tid)
lIdx = [];
if isempty(TLine) || height(TLine) == 0, return; end
tok = regexp(tid, '^(-?\d+)\s*-\s*(-?\d+)$', 'tokens', 'once');
if ~isempty(tok)
    fb = str2double(tok{1});
    tb = str2double(tok{2});
    lIdx = find((TLine.fbus == fb & TLine.tbus == tb) | (TLine.fbus == tb & TLine.tbus == fb), 1);
    return;
end
num = str2double(regexprep(tid, '^#', ''));
if ~isnan(num) && num >= 1 && num <= height(TLine)
    lIdx = round(num);
end
end

function rCh = local_ref_row(ch)
if contains(ch, 'dw') || contains(ch, 'omega') || contains(ch, 'freq')
    rCh = 1;
elseif contains(ch, 'dv') || contains(ch, 'volt')
    rCh = 3;
elseif contains(ch, 'dq')
    rCh = 4;
else
    rCh = 2; % dP_ref_pu default
end
end
