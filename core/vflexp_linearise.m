function linRes = vflexp_linearise(session, forceReinit)
%VFLEXP_LINEARISE Linearise VFlexP around the initialized steady-state operating point.
%
%   linRes = vflexp_linearise(session, forceReinit)
%
%   Uses in-memory linio objects (NEVER setlinio, keeping VFlexP.slx
%   Dirty == 'off') to extract the full state-space model (A, B, C, D)
%   across AC/DC network and agent perturbation/measurement channels, and
%   maps Simulink block paths to human-readable labels via
%   vflexp_map_state_names.

if nargin < 2 || isempty(forceReinit)
    forceReinit = false;
end

tStart = tic;
if forceReinit || isempty(session.initResult)
    vflexp_build_and_init(session);
end

mdlWks   = get_param('VFlexP', 'ModelWorkspace');
Nbuses   = mdlWks.getVariable('Nbuses');
buses_lf = mdlWks.getVariable('buses_lf');
Nag      = numel(buses_lf);
DC       = mdlWks.getVariable('DC');
hasDC    = isstruct(DC) && isfield(DC, 'include') && logical(DC.include) && DC.Nbuses > 0;

% Ensure disturbance time-series are zeroed during t=0 linearisation
savedEvents = session.eventList;
if ~isempty(savedEvents) && height(savedEvents) > 0
    session.eventList = VFlexPSession.emptyEventsTable();
    evZero = vflexp_build_events_ws(session);
    local_assign_ev_ws(mdlWks, evZero);
    cEv = onCleanup(@() local_restore_events(session, mdlWks, savedEvents));
    if isempty(cEv), return; end
end

% Build in-memory linearization I/O points (zero model mutation)
io = [ ...
    linio('VFlexP/Grid/in_current/u_bus_P', 1, 'input'); ...
    linio('VFlexP/Grid/in_current/u_bus_Q', 1, 'input')];
if hasDC
    io(end+1, 1) = linio('VFlexP/Grid/u_dc_I', 1, 'input');
end
if Nag > 0
    io(end+1, 1) = linio('VFlexP/Grid/u_gen_ref', 1, 'input');
end

io(end+1, 1) = linio('VFlexP/Grid/y_vac', 1, 'output');
if hasDC
    io(end+1, 1) = linio('VFlexP/Grid/y_vdc', 1, 'output');
end
if Nag > 0
    io(end+1, 1) = linio('VFlexP/Grid/y_elec', 1, 'output');
end

opt = linearizeOptions( ...
    'UseFullBlockNameLabels', 'on', ...
    'SampleTime', 0);

linsys = linearize('VFlexP', io, opt);

% Map input channel column indices cleanly
inMap = struct();
colCursor = 0;
inMap.bus_P = (colCursor + 1):(colCursor + Nbuses);
colCursor   = colCursor + Nbuses;
inMap.bus_Q = (colCursor + 1):(colCursor + Nbuses);
colCursor   = colCursor + Nbuses;

if hasDC
    inMap.dc_I = (colCursor + 1):(colCursor + DC.Nbuses);
    colCursor  = colCursor + DC.Nbuses;
else
    inMap.dc_I = zeros(1, 0);
end

if Nag > 0
    % u_gen_ref is [6 x Nag] serialized column-major: rows 1..4 = [omega*, P*, V*, Q*]
    baseIdx = colCursor + (0:(Nag - 1)) * 6;
    inMap.gen_omega = baseIdx + 1;
    inMap.gen_P     = baseIdx + 2;
    inMap.gen_V     = baseIdx + 3;
    inMap.gen_Q     = baseIdx + 4;
else
    inMap.gen_omega = zeros(1, 0);
    inMap.gen_P     = zeros(1, 0);
    inMap.gen_V     = zeros(1, 0);
    inMap.gen_Q     = zeros(1, 0);
end

% Map output channel row indices cleanly
outMap = struct();
rowCursor = 0;
% y_vac is [3 x Nbuses]; keep rows 1 & 2 (vd, vq), skip inert zero-sequence row 3
vacBase = rowCursor + (0:(Nbuses - 1)) * 3;
outMap.vac = sort([vacBase + 1, vacBase + 2]);
rowCursor  = rowCursor + 3 * Nbuses;

if hasDC
    outMap.vdc = (rowCursor + 1):(rowCursor + DC.Nbuses);
    rowCursor  = rowCursor + DC.Nbuses;
else
    outMap.vdc = zeros(1, 0);
end

if Nag > 0
    % y_elec is [6 x Nag]; rows 1..4 = [omega, P, V, Q]
    elecBase = rowCursor + (0:(Nag - 1)) * 6;
    outMap.elec   = sort([elecBase + 1, elecBase + 2, elecBase + 3, elecBase + 4]);
    outMap.elec_f = elecBase + 1;
    outMap.elec_P = elecBase + 2;
    outMap.elec_V = elecBase + 3;
    outMap.elec_Q = elecBase + 4;
else
    outMap.elec   = zeros(1, 0);
    outMap.elec_f = zeros(1, 0);
    outMap.elec_P = zeros(1, 0);
    outMap.elec_V = zeros(1, 0);
    outMap.elec_Q = zeros(1, 0);
end

[stateLabels, stateMeta] = vflexp_map_state_names(session, linsys.StateName);

% Truncate invariant reference-frame angle state(s) (d(theta_ref)/dt == 0)
elimIdx = local_find_ref_angle_states(session, mdlWks, linsys, stateMeta);
if ~isempty(elimIdx)
    keepIdx = setdiff((1:size(linsys.A, 1))', elimIdx(:), 'stable');
    linsys  = modred(linsys, elimIdx, 'Truncate');
    stateLabels = stateLabels(keepIdx);
    stateMeta   = stateMeta(keepIdx, :);
    stateMeta.Index = (1:numel(keepIdx))';
end

elapsed = toc(tStart);
linRes = struct( ...
    'linsys', linsys, ...
    'A', double(linsys.A), ...
    'B', double(linsys.B), ...
    'C', double(linsys.C), ...
    'D', double(linsys.D), ...
    'stateLabels', stateLabels, ...
    'stateMeta', stateMeta, ...
    'inMap', inMap, ...
    'outMap', outMap, ...
    'nStates', size(linsys.A, 1), ...
    'elapsed', elapsed);

session.linResult = linRes;
session.log('Linearised VFlexP in %.3f s (%d states, %d inputs, %d outputs).', ...
    elapsed, size(linsys.A, 1), size(linsys.B, 2), size(linsys.C, 1));
end

function elimIdx = local_find_ref_angle_states(session, mdlWks, linsys, stateMeta)
elimIdx = zeros(0, 1);
nX = size(linsys.A, 1);
if nX == 0 || isempty(stateMeta)
    return;
end
try
    fref = mdlWks.getVariable('fref');
    useSync = mdlWks.getVariable('use_sync_frame');
catch
    return;
end
if ~isstruct(fref) || ~isfield(fref, 'ref_units')
    return;
end
nAreas = numel(fref.ref_units);
if isscalar(useSync)
    useSync = repmat(double(useSync), 1, nAreas);
else
    useSync = reshape(double(useSync), 1, []);
end

[~, TgenAC] = session.splitDomains();
if isempty(TgenAC) || height(TgenAC) == 0
    return;
end

A = double(linsys.A);
B = double(linsys.B);
for a = 1:min(nAreas, numel(useSync))
    if useSync(a) >= 0.5
        continue;
    end
    gRow = round(fref.ref_units(a));
    if gRow < 1 || gRow > height(TgenAC)
        continue;
    end
    refBus  = sprintf('%d', round(TgenAC.bus(gRow)));
    refType = string(TgenAC.type(gRow));
    cand = find(stateMeta.ElementKind == "ac_gen" & ...
                stateMeta.UnitType == refType & ...
                stateMeta.ElementId == string(refBus) & ...
                (endsWith(stateMeta.RawPath, '/theta') | ...
                 contains(stateMeta.RawPath, '/theta(') | ...
                 endsWith(stateMeta.RawPath, '/frame') | ...
                 contains(stateMeta.RawPath, '/frame(')));
    for c = cand(:)'
        if norm(A(c, :), inf) < 1e-10 && (isempty(B) || norm(B(c, :), inf) < 1e-10)
            elimIdx(end+1, 1) = c; %#ok<AGROW>
        end
    end
end
elimIdx = unique(elimIdx, 'stable');
end

function local_assign_ev_ws(mdlWks, evWs)
fn = fieldnames(evWs);
for k = 1:numel(fn)
    assignin(mdlWks, fn{k}, evWs.(fn{k}));
end
end

function local_restore_events(session, mdlWks, savedEvents)
session.eventList = savedEvents;
evWs = vflexp_build_events_ws(session);
local_assign_ev_ws(mdlWks, evWs);
end
