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

nGenAC = height(TgenAC);
nGenDC = height(TgenDC);

for r = 1:height(sigTbl)
    t0  = max(0.0, double(sigTbl.Time(r)));
    dur = double(sigTbl.Duration(r));
    mag = double(sigTbl.Magnitude(r));
    tt  = lower(strtrim(char(string(sigTbl.TargetType(r)))));
    ch  = lower(strtrim(char(string(sigTbl.Channel(r)))));
    tid = strtrim(char(string(sigTbl.TargetId(r))));

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

evWs = struct();
evWs.ev_bus_dP   = struct('time', tVec, 'signals', struct('values', valBusP, 'dimensions', [1, nAC]));
evWs.ev_bus_dQ   = struct('time', tVec, 'signals', struct('values', valBusQ, 'dimensions', [1, nAC]));
evWs.ev_dc_dI    = struct('time', tVec, 'signals', struct('values', valDcI,  'dimensions', [1, nDC]));
evWs.ev_gen_dref = struct('time', tVec, 'signals', struct('values', valRef,  'dimensions', [6, nAg]));
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
evWs.t_close_lines = t_close_lines;
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
