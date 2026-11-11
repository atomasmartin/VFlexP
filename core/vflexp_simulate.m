function simRes = vflexp_simulate(session, reinit)
%VFLEXP_SIMULATE Run time-domain simulation of VFlexP with zero .slx mutation.
%
%   simRes = vflexp_simulate(session)
%   simRes = vflexp_simulate(session, reinit)

if nargin < 2 || isempty(reinit)
    reinit = isempty(session.initResult);
end

if reinit
    vflexp_build_and_init(session);
else
    % Refresh event FromWorkspace structs in ModelWorkspace in case events or t_end changed
    mdlWks = get_param('VFlexP', 'ModelWorkspace');
    evWs = vflexp_build_events_ws(session);
    assignin(mdlWks, 'ev_bus_dP', evWs.ev_bus_dP);
    assignin(mdlWks, 'ev_bus_dQ', evWs.ev_bus_dQ);
    assignin(mdlWks, 'ev_dc_dI', evWs.ev_dc_dI);
    assignin(mdlWks, 'ev_gen_dref', evWs.ev_gen_dref);
    assignin(mdlWks, 't_close_lines', evWs.t_close_lines);
    if any(evWs.t_close_lines > 0)
        assignin(mdlWks, 'ac_grid_mode', 1);
    end
end

tStart = tic;
scfg = session.config.sim;
inSim = Simulink.SimulationInput('VFlexP');
inSim = inSim.setModelParameter( ...
    'StartTime', '0', ...
    'StopTime', num2str(double(scfg.t_end), '%.12g'), ...
    'Solver', char(string(scfg.method)), ...
    'RelTol', num2str(double(scfg.rtol), '%.12g'), ...
    'AbsTol', num2str(double(scfg.atol), '%.12g'));
if isnumeric(scfg.max_step) && ~isempty(scfg.max_step) && isfinite(scfg.max_step) && scfg.max_step > 0
    inSim = inSim.setModelParameter('MaxStep', num2str(double(scfg.max_step), '%.12g'));
else
    inSim = inSim.setModelParameter('MaxStep', 'auto');
end

simOut = sim(inSim);
elapsed = toc(tStart);

[TbusAC, TgenAC, ~, TbusDC, TgenDC, ~, Tintconv, TFb, Sb] = session.splitDomains();
mdlWks = get_param('VFlexP', 'ModelWorkspace');
fref   = mdlWks.getVariable('fref');

% Extract logged signals (timeseries or struct-with-time)
[tVac, vacData]   = local_extract_3d(simOut, 'vflexp_out_vac');   % [3 x Nbuses x Nt]
[~,    vdcData]   = local_extract_2d(simOut, 'vflexp_out_vdc');   % [Nt x Ndc]
[tElec, elecData] = local_extract_3d(simOut, 'vflexp_out_elec');  % [6 x Nag x Nt]

if ~isempty(tVac)
    t = tVac(:);
else
    t = tElec(:);
end
nT = numel(t);

vd = permute(vacData(1, :, :), [3, 2, 1]); % [Nt x Nbuses]
vq = permute(vacData(2, :, :), [3, 2, 1]);
vac_mag     = hypot(vd, vq);
vac_ang_deg = rad2deg(atan2(vq, vd));

omega_pu = permute(elecData(1, :, :), [3, 2, 1]); % [Nt x Nag]
P_pu     = permute(elecData(2, :, :), [3, 2, 1]);
Vgen_pu  = permute(elecData(3, :, :), [3, 2, 1]);
Q_pu     = permute(elecData(4, :, :), [3, 2, 1]);

fbAreas = double(TFb.Fb(:))';
fbAgents = fbAreas(fref.Fb_units);
omega_Hz = omega_pu .* repmat(fbAgents, nT, 1);
P_MW     = P_pu * Sb;
Q_MVAr   = Q_pu * Sb;

if isempty(TbusDC) || height(TbusDC) == 0
    vdc = zeros(nT, 0);
    dcBusLabels = strings(0, 1);
else
    vdc = vdcData;
    dcBusLabels = compose("DC Bus %d", TbusDC.bus_i(:));
end

busLabels = compose("Bus %d", TbusAC.bus_i(:));
agentLabels = strings(0, 1);
for g = 1:height(TgenAC)
    agentLabels(end+1, 1) = sprintf('%s @ Bus %d', TgenAC.type(g), TgenAC.bus(g)); %#ok<AGROW>
end
for g = 1:height(TgenDC)
    agentLabels(end+1, 1) = sprintf('%s @ DC Bus %d', TgenDC.type(g), TgenDC.bus(g)); %#ok<AGROW>
end
for c = 1:height(Tintconv)
    agentLabels(end+1, 1) = sprintf('%s (%d-%d)', Tintconv.type(c), Tintconv.fbus(c), Tintconv.tbus(c)); %#ok<AGROW>
end

simRes = struct( ...
    't', t, ...
    'vac_mag', vac_mag, ...
    'vac_ang_deg', vac_ang_deg, ...
    'vac_d', vd, ...
    'vac_q', vq, ...
    'vdc', vdc, ...
    'omega_pu', omega_pu, ...
    'omega_Hz', omega_Hz, ...
    'P_MW', P_MW, ...
    'Q_MVAr', Q_MVAr, ...
    'Vgen_pu', Vgen_pu, ...
    'busLabels', busLabels, ...
    'dcBusLabels', dcBusLabels, ...
    'agentLabels', agentLabels, ...
    'elapsed', elapsed);

session.simResult = simRes;
session.log('Simulated [0, %.2f] s in %.3f s (%d time steps, solver = %s).', ...
    double(scfg.t_end), elapsed, nT, char(string(scfg.method)));
end

% -------------------------------------------------------------------------
function [t, arr3] = local_extract_3d(simOut, varName)
raw = simOut.get(varName);
if isa(raw, 'timeseries')
    t = raw.Time(:);
    d = double(raw.Data);
    nT = numel(t);
    if ndims(d) == 3
        if size(d, 3) == nT
            arr3 = d;
        else
            arr3 = permute(d, [2, 3, 1]);
        end
    else
        % 2D when 2nd dim is 1: Simulink may store [Nt x rows] or [rows x Nt]
        if size(d, 1) == nT
            arr3 = reshape(d', size(d, 2), 1, nT);
        else
            arr3 = reshape(d, size(d, 1), 1, nT);
        end
    end
elseif isstruct(raw)
    t = raw.time(:);
    v = double(raw.signals.values);
    nT = numel(t);
    if ndims(v) == 3
        arr3 = v;
    elseif size(v, 1) == nT
        arr3 = reshape(v', size(v, 2), 1, nT);
    else
        arr3 = reshape(v, size(v, 1), 1, nT);
    end
else
    t = simOut.tout(:);
    arr3 = double(raw);
end
end

function [t, mat2] = local_extract_2d(simOut, varName)
raw = simOut.get(varName);
if isa(raw, 'timeseries')
    t = raw.Time(:);
    d = squeeze(double(raw.Data));
    nT = numel(t);
    if isvector(d)
        mat2 = d(:);
    elseif size(d, 1) == nT
        mat2 = d;
    else
        mat2 = d';
    end
elseif isstruct(raw)
    t = raw.time(:);
    v = squeeze(double(raw.signals.values));
    nT = numel(t);
    if isvector(v)
        mat2 = v(:);
    elseif size(v, 1) == nT
        mat2 = v;
    else
        mat2 = v';
    end
else
    t = simOut.tout(:);
    mat2 = double(raw);
end
end
