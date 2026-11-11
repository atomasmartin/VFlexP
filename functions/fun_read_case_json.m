function settings = fun_read_case_json(json_file, gcb)
%FUN_READ_CASE_JSON Read a VFlexP/VFlexPy <case>.json configuration file.
%
%   settings = fun_read_case_json(json_file)
%   settings = fun_read_case_json(json_file, gcb)
%
%   Reads simulation and power-flow configuration from a JSON file
%   (t_end, solver.rtol, solver.atol, solver.max_step, powerflow.*) and,
%   when gcb is provided, applies StopTime / RelTol / AbsTol to the
%   active Simulink model and stores the struct in the ModelWorkspace.

settings = struct( ...
    't_end', [], ...
    'frame_policy', 'first-gen', ...
    'events', {strings(0, 1)}, ...
    'solver', struct('method', 'ode15s', 'rtol', 1e-6, 'atol', 1e-8, 'max_step', 'auto'), ...
    'powerflow', struct('tolerance', 1e-8, 'max_iterations', 50, 'enforce_q_limits', false));

if nargin < 1 || isempty(json_file) || ~isfile(json_file)
    return
end

raw = jsondecode(fileread(json_file));

if isfield(raw, 'model') && isstruct(raw.model)
    if isfield(raw.model, 'frame_policy') && ~isempty(raw.model.frame_policy)
        settings.frame_policy = char(strtrim(string(raw.model.frame_policy)));
    end
elseif isfield(raw, 'frame_policy') && ~isempty(raw.frame_policy)
    settings.frame_policy = char(strtrim(string(raw.frame_policy)));
end

% Support both flat keys (t_end, solver, powerflow) and VFlexPy action-keyed
% blocks ("dyn.simulate" -> dyn_simulate, "pf.run" -> pf_run)
hasExplicitRtol = false;
hasExplicitAtol = false;
if isfield(raw, 'dyn_simulate') && isstruct(raw.dyn_simulate)
    ds = raw.dyn_simulate;
    if isfield(ds, 't_end') && ~isempty(ds.t_end)
        settings.t_end = double(ds.t_end);
    end
    if isfield(ds, 'method') && ~isempty(ds.method)
        settings.solver.method = local_map_solver(char(string(ds.method)));
    end
    if isfield(ds, 'rtol') && ~isempty(ds.rtol)
        settings.solver.rtol = double(ds.rtol);
        hasExplicitRtol = true;
    end
    if isfield(ds, 'atol') && ~isempty(ds.atol)
        settings.solver.atol = double(ds.atol);
        hasExplicitAtol = true;
    end
    if isfield(ds, 'max_step') && ~isempty(ds.max_step)
        settings.solver.max_step = ds.max_step;
    elseif isfield(ds, 'step') && ~isempty(ds.step)
        settings.solver.max_step = ds.step;
    end
    if isfield(ds, 'event') && ~isempty(ds.event)
        settings.events = reshape(string(ds.event), [], 1);
    end
    if isfield(ds, 'frame_policy') && ~isempty(ds.frame_policy)
        settings.frame_policy = char(strtrim(string(ds.frame_policy)));
    end
end

if isfield(raw, 't_end') && ~isempty(raw.t_end)
    settings.t_end = double(raw.t_end);
end
if isfield(raw, 'event') && ~isempty(raw.event)
    settings.events = reshape(string(raw.event), [], 1);
end

if isfield(raw, 'solver') && isstruct(raw.solver)
    s = raw.solver;
    if isfield(s, 'method') && ~isempty(s.method)
        settings.solver.method = local_map_solver(char(string(s.method)));
    end
    if isfield(s, 'rtol') && ~isempty(s.rtol)
        settings.solver.rtol = double(s.rtol);
        hasExplicitRtol = true;
    end
    if isfield(s, 'atol') && ~isempty(s.atol)
        settings.solver.atol = double(s.atol);
        hasExplicitAtol = true;
    end
    if isfield(s, 'max_step') && ~isempty(s.max_step)
        settings.solver.max_step = s.max_step;
    end
end

if hasExplicitRtol && ~hasExplicitAtol && settings.solver.rtol > 1e-5
    settings.solver.atol = max(1e-8, 1e-3 * settings.solver.rtol);
end

for pf_field = {'pf_run', 'powerflow'}
    fn = pf_field{1};
    if isfield(raw, fn) && isstruct(raw.(fn))
        pf = raw.(fn);
        if isfield(pf, 'solver') && ~isempty(pf.solver)
            settings.powerflow.solver = char(string(pf.solver));
        end
        if isfield(pf, 'tolerance') && ~isempty(pf.tolerance)
            settings.powerflow.tolerance = double(pf.tolerance);
        end
        if isfield(pf, 'max_iterations') && ~isempty(pf.max_iterations)
            settings.powerflow.max_iterations = double(pf.max_iterations);
        end
        if isfield(pf, 'enforce_q_limits') && ~isempty(pf.enforce_q_limits)
            settings.powerflow.enforce_q_limits = logical(pf.enforce_q_limits);
        end
    end
end

if nargin >= 2 && ~isempty(gcb)
    try
        mdl = bdroot(gcb);
        if ~bdIsLibrary(mdl)
            if ~isempty(settings.t_end)
                set_param(mdl, 'StopTime', num2str(settings.t_end, '%.12g'));
            end
            has_solver_cfg = (isfield(raw, 'solver') && isstruct(raw.solver)) || ...
                             (isfield(raw, 'dyn_simulate') && isstruct(raw.dyn_simulate));
            if has_solver_cfg
                if ~isempty(settings.solver.rtol)
                    set_param(mdl, 'RelTol', num2str(settings.solver.rtol, '%.12g'));
                end
                if ~isempty(settings.solver.atol)
                    set_param(mdl, 'AbsTol', num2str(settings.solver.atol, '%.12g'));
                end
                if isnumeric(settings.solver.max_step) && ~isempty(settings.solver.max_step)
                    set_param(mdl, 'MaxStep', num2str(double(settings.solver.max_step), '%.12g'));
                end
            end
            mdlWks = get_param(mdl, 'ModelWorkspace');
            bvarname = get_param(gcb, 'bvarname');
            assignin(mdlWks, [bvarname 'case_settings'], settings);
        end
    catch
    end
end
end

function mOut = local_map_solver(mIn)
mClean = lower(strtrim(char(string(mIn))));
switch mClean
    case {'radau', 'bdf2'}
        mOut = 'ode23tb';
    case {'bdf', 'lsoda', 'auto'}
        mOut = 'ode15s';
    case 'rk45'
        mOut = 'ode45';
    case 'dop853'
        mOut = 'ode89';
    otherwise
        mOut = strtrim(char(string(mIn)));
end
end
