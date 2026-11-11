function initRes = vflexp_build_and_init(session)
%VFLEXP_BUILD_AND_INIT Populate VFlexP ModelWorkspace and initialize all grid & unit states.
%
%   initRes = vflexp_build_and_init(session)
%
%   Writes exclusively to ModelWorkspace of VFlexP.slx (zero set_param on
%   VFlexP blocks, preserving LinkStatus == 'resolved' and Dirty == 'off').
%   Uses exact analytical steady-state initialization for AC/DC grids, GFM,
%   intconv_GFL, and intconv_VSM, and an in-memory temporary Simulink
%   harness (findop) for synchronous machines and remaining units.

tStart = tic;
if isempty(session.pfResult) || ~session.pfResult.converged
    vflexp_run_pf(session);
    if ~session.pfResult.converged
        error('VFlexP:BuildAndInit:PFNotConverged', ...
            'Power Flow must converge before initializing the dynamic model.');
    end
end

session.evalTables();
[TbusAC, TgenAC, TlineAC, TbusDC, TgenDC, TlineDC, Tintconv, TFb, Sb] = session.splitDomains();

if isempty(TbusAC) || height(TbusAC) == 0
    error('VFlexP:BuildAndInit:NoACBuses', 'Case must contain at least one AC bus.');
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
slxPath  = fullfile(repoRoot, 'Simulink', 'VFlexP.slx');
if bdIsLoaded('VFlexP')
    if ~strcmpi(get_param('VFlexP', 'FileName'), slxPath)
        close_system('VFlexP', 0);
    end
end
if ~bdIsLoaded('VFlexP')
    load_system(slxPath);
end
mdlWks = get_param('VFlexP', 'ModelWorkspace');
manifest = local_get_manifest(repoRoot);

% Reset clean workspace defaults
evalin(mdlWks, 'vflexp_ws_defaults;');

% -------------------------------------------------------------------------
% 1. Base frequencies & area reference mapping (fref)
% -------------------------------------------------------------------------
if isempty(TFb) || height(TFb) == 0
    TFb = table(1, 50, 'VariableNames', {'Area', 'Fb'});
end
omega_base = 2 * pi * reshape(double(TFb.Fb), 1, []);
assignin(mdlWks, 'Sb', Sb);
assignin(mdlWks, 'Fb', double(TFb.Fb(1)));
assignin(mdlWks, 'omega_base', omega_base);

fref = struct();
fref.areas = reshape(double(TFb.Area), 1, []);
nAreas = numel(fref.areas);
[fref.ref_units, use_sync_frame] = local_resolve_area_frames(fref.areas, TbusAC, TgenAC, session);
assignin(mdlWks, 'use_sync_frame', use_sync_frame);

Nbuses = height(TbusAC);
fref.buses = ones(1, Nbuses);
for b = 1:Nbuses
    pos = find(fref.areas == TbusAC.area(b), 1);
    if ~isempty(pos), fref.buses(b) = pos; end
end

hasACLines = ~isempty(TlineAC) && height(TlineAC) > 0;
if hasACLines
    nLinesAC = height(TlineAC);
    fref.lines = ones(1, nLinesAC);
    for l = 1:nLinesAC
        bPos = find(TbusAC.bus_i == TlineAC.fbus(l), 1);
        if ~isempty(bPos)
            aPos = find(fref.areas == TbusAC.area(bPos), 1);
            if ~isempty(aPos), fref.lines(l) = aPos; end
        end
    end
else
    fref.lines = 1;
end

loadMaskAC = (TbusAC.Pd ~= 0) | (TbusAC.Qd ~= 0);
loadIdxAC  = find(loadMaskAC);
if ~isempty(loadIdxAC)
    fref.loads = fref.buses(loadIdxAC);
else
    fref.loads = 1;
end

% -------------------------------------------------------------------------
% 2. AC network parameters, complexity vectors & operating point
% -------------------------------------------------------------------------
bus_ids = TbusAC.bus_i(:);
v_buses_mod   = reshape(double(TbusAC.Vm), 1, []);
v_buses_angle = reshape(deg2rad(double(TbusAC.Va)), 1, []);
vd_buses = real(v_buses_mod .* exp(1i * v_buses_angle));
vq_buses = imag(v_buses_mod .* exp(1i * v_buses_angle));

assignin(mdlWks, 'Nbuses', Nbuses);
assignin(mdlWks, 'bus', bus_ids);
assignin(mdlWks, 'v_buses_mod', v_buses_mod);

gs = reshape(double(TbusAC.Gs), 1, []);
gsSafe = gs;
gsSafe(gsSafe == 0) = 1e-12;
Rcf = Sb ./ gsSafe;
assignin(mdlWks, 'Rcf', Rcf);

C_to_pu = reshape(double(TbusAC.Bs), 1, []) / Sb;

if hasACLines
    nL = height(TlineAC);
    from_idx = zeros(1, nL);
    to_idx   = zeros(1, nL);
    for l = 1:nL
        from_idx(l) = find(bus_ids == TlineAC.fbus(l), 1);
        to_idx(l)   = find(bus_ids == TlineAC.tbus(l), 1);
    end
    R_Lines = reshape(double(TlineAC.r), 1, []);
    L_Lines = reshape(double(TlineAC.x), 1, []);
    C_Lines = reshape(double(TlineAC.b), 1, []);
    tap_lines = reshape(double(TlineAC.ratio), 1, []);
    tap_lines(tap_lines == 0) = 1.0;
    shift_lines = reshape(deg2rad(double(TlineAC.angle)), 1, []);

    Lines_fromM = zeros(Nbuses, nL);
    Lines_toM   = zeros(Nbuses, nL);
    for l = 1:nL
        Lines_fromM(from_idx(l), l) = 1;
        Lines_toM(to_idx(l), l)     = 1;
    end

    Cbuses = C_to_pu;
    for b = 1:Nbuses
        Cbuses(b) = C_to_pu(b) + 0.5 * sum(C_Lines(from_idx == b | to_idx == b));
    end

    line_complexity = local_encode_ac_line_complexity(TlineAC.Complexity);
    minTau = 1e-5;
    for l = 1:nL
        if line_complexity(l) == 1
            wbL = omega_base(fref.lines(l));
            if L_Lines(l) <= 0 || (R_Lines(l) > 0 && (L_Lines(l) / (R_Lines(l) * wbL)) < minTau)
                line_complexity(l) = 3;
            end
        end
    end

    vfrom = v_buses_mod(from_idx) .* exp(1i * v_buses_angle(from_idx));
    vfrom = (vfrom ./ tap_lines) .* exp(-1i * shift_lines);
    vto   = v_buses_mod(to_idx)   .* exp(1i * v_buses_angle(to_idx));
    dv_d  = real(vfrom - vto);
    dv_q  = imag(vfrom - vto);
    z2    = max(R_Lines.^2 + L_Lines.^2, 1e-12);
    i_lines_d = (R_Lines .* dv_d + L_Lines .* dv_q) ./ z2;
    i_lines_q = (R_Lines .* dv_q - L_Lines .* dv_d) ./ z2;
else
    from_idx    = 1;
    to_idx      = 1;
    R_Lines     = 1;
    L_Lines     = 1;
    C_Lines     = 0;
    tap_lines   = 1;
    shift_lines = 0;
    Lines_fromM = zeros(Nbuses, 1);
    Lines_toM   = zeros(Nbuses, 1);
    Cbuses      = C_to_pu;
    line_complexity = 0;
    i_lines_d   = 0;
    i_lines_q   = 0;
end

bus_complexity = local_encode_ac_bus_complexity(TbusAC.Complexity_bus);
for b = 1:Nbuses
    if bus_complexity(b) == 1
        wbB = omega_base(fref.buses(b));
        gBus = gs(b) / Sb;
        if Cbuses(b) <= 0 || (gBus > 0 && (Cbuses(b) / (gBus * wbB)) < 1e-5)
            bus_complexity(b) = 3;
        end
    end
end
Cbuses(bus_complexity == 1 & Cbuses <= 0) = 1e-12;

assignin(mdlWks, 'From_buses', from_idx);
assignin(mdlWks, 'To_buses', to_idx);
assignin(mdlWks, 'R_Lines', R_Lines);
assignin(mdlWks, 'L_Lines', L_Lines);
assignin(mdlWks, 'C_Lines', C_Lines);
assignin(mdlWks, 'tap_lines', tap_lines);
assignin(mdlWks, 'shift_lines', shift_lines);
assignin(mdlWks, 'Lines_fromM', sparse(Lines_fromM));
assignin(mdlWks, 'Lines_toM', sparse(Lines_toM));
assignin(mdlWks, 'Cbuses', Cbuses);
assignin(mdlWks, 'bus_complexity', bus_complexity);
assignin(mdlWks, 'line_complexity', line_complexity);
for kComp = 1:4
    assignin(mdlWks, sprintf('has_bus_c%d', kComp), double(any(bus_complexity == kComp)));
end
for kComp = 1:3
    assignin(mdlWks, sprintf('has_line_c%d', kComp), double(hasACLines && any(line_complexity == kComp)));
end

if ~isempty(loadIdxAC)
    Ploads = reshape(double(TbusAC.Pd(loadIdxAC)), 1, []);
    Qloads = reshape(double(TbusAC.Qd(loadIdxAC)), 1, []);
    vmL    = v_buses_mod(loadIdxAC);
    s2L    = max(Ploads.^2 + Qloads.^2, 1e-20);
    R_Loads = (vmL.^2) .* Sb .* (Ploads ./ s2L);
    L_Loads = (vmL.^2) .* Sb .* (Qloads ./ s2L);
    % Guard pure P or pure Q loads when series/parallel RL blocks evaluate
    R_Loads(R_Loads == 0 & L_Loads == 0) = 1e6;
    load_complexity = local_encode_ac_load_complexity(TbusAC.Complexity_load(loadIdxAC));
    for kLd = 1:numel(load_complexity)
        bIdxLd = loadIdxAC(kLd);
        wbLd   = omega_base(fref.loads(kLd));
        if load_complexity(kLd) == 1
            tauFast = (L_Loads(kLd) <= 0) || ...
                      (R_Loads(kLd) > 0 && (L_Loads(kLd) / (R_Loads(kLd) * wbLd)) < 1e-5);
            if tauFast
                load_complexity(kLd) = 3;
            end
        elseif load_complexity(kLd) == 6
            tauFast = (Qloads(kLd) <= 0) || ...
                      (Ploads(kLd) > 0 && ((Qloads(kLd) / Ploads(kLd)) / wbLd) < 1e-5);
            if tauFast
                if bus_complexity(bIdxLd) == 1
                    load_complexity(kLd) = 3;
                else
                    Qloads(kLd) = 1e-3 * Sb;
                end
            end
        end
    end

    vdL = vd_buses(loadIdxAC);
    vqL = vq_buses(loadIdxAC);
    z2Ld = max(R_Loads.^2 + L_Loads.^2, 1e-12);
    i_loads_d = (R_Loads .* vdL + L_Loads .* vqL) ./ z2Ld;
    i_loads_q = (R_Loads .* vqL - L_Loads .* vdL) ./ z2Ld;
    Load_buses = loadIdxAC(:);
else
    Load_buses = 1;
    Ploads     = 0;
    Qloads     = 0;
    R_Loads    = 1e6;
    L_Loads    = 1e6;
    load_complexity = 0;
    i_loads_d  = 0;
    i_loads_q  = 0;
end
assignin(mdlWks, 'Load_buses', Load_buses);
assignin(mdlWks, 'Ploads', Ploads);
assignin(mdlWks, 'Qloads', Qloads);
assignin(mdlWks, 'R_Loads', R_Loads);
assignin(mdlWks, 'L_Loads', L_Loads);
assignin(mdlWks, 'load_complexity', load_complexity);
for kComp = 1:9
    assignin(mdlWks, sprintf('has_load_c%d', kComp), double(~isempty(loadIdxAC) && any(load_complexity == kComp)));
end

% Admittance matrices for Static AC grid mode
[Y_ac, Y_ac_noloads] = local_admitt(Sb, TbusAC, TlineAC);
assignin(mdlWks, 'admitt_matrix', Y_ac);
assignin(mdlWks, 'iadmitt_matrix', local_safe_inv(Y_ac));
assignin(mdlWks, 'admitt_matrix_noloads', Y_ac_noloads);
assignin(mdlWks, 'iadmitt_matrix_noloads', local_safe_inv(Y_ac_noloads));

% Events workspace variables
evWs = vflexp_build_events_ws(session);
assignin(mdlWks, 'ev_bus_dP', evWs.ev_bus_dP);
assignin(mdlWks, 'ev_bus_dQ', evWs.ev_bus_dQ);
assignin(mdlWks, 'ev_dc_dI', evWs.ev_dc_dI);
assignin(mdlWks, 'ev_gen_dref', evWs.ev_gen_dref);
assignin(mdlWks, 'ev_bus_constZ', evWs.ev_bus_constZ);
assignin(mdlWks, 't_close_lines', evWs.t_close_lines);

allStaticAC = all(bus_complexity == 3) && ...
              all(load_complexity == 3 | load_complexity == 0) && ...
              all(line_complexity == 3 | line_complexity == 0) && ...
              all(evWs.t_close_lines == 0);
ac_grid_mode = double(~allStaticAC);
assignin(mdlWks, 'ac_grid_mode', ac_grid_mode);

% Operating point for Dynamic AC Grid
opGridAC = mdlWks.getVariable('op_Grid_AC_Dynamic');
opGridAC.Buses_dynamic_Active_Impl_vd = local_nonempty(vd_buses(bus_complexity == 1));
opGridAC.Buses_dynamic_Active_Impl_vq = local_nonempty(vq_buses(bus_complexity == 1));
if any(bus_complexity == 3)
    nB3 = sum(bus_complexity == 3);
    icB3 = [vd_buses(bus_complexity == 3); vq_buses(bus_complexity == 3); zeros(1, nB3)];
    opGridAC.Buses_static_omega_Active_Impl_LPF_I = icB3;
    opGridAC.bus_static_omega_ic = icB3;
end
if hasACLines && any(line_complexity == 1)
    opGridAC.Lines_dynamic_Active_Impl_i_d = i_lines_d(line_complexity == 1);
    opGridAC.Lines_dynamic_Active_Impl_i_q = i_lines_q(line_complexity == 1);
end
if ~isempty(loadIdxAC)
    if any(load_complexity == 1)
        opGridAC.Loads_series_RL_dynamic_Active_Impl_i_d = i_loads_d(load_complexity == 1);
        opGridAC.Loads_series_RL_dynamic_Active_Impl_i_q = i_loads_q(load_complexity == 1);
    end
    if any(load_complexity == 5)
        opGridAC.Loads_dynamic_PQ_Active_Impl_i_d_I = i_loads_d(load_complexity == 5);
        opGridAC.Loads_dynamic_PQ_Active_Impl_i_q_I = i_loads_q(load_complexity == 5);
    end
    for parCode = [6, 7, 8]
        if any(load_complexity == parCode)
            lbPar = loadIdxAC(load_complexity == parCode);
            qPar  = Qloads(load_complexity == parCode);
            idPar =  vq_buses(lbPar) .* qPar ./ (v_buses_mod(lbPar).^2) / Sb;
            iqPar = -vd_buses(lbPar) .* qPar ./ (v_buses_mod(lbPar).^2) / Sb;
            switch parCode
                case 6
                    opGridAC.Loads_parallel_RL_dynamic_Active_Impl_i_d = idPar;
                    opGridAC.Loads_parallel_RL_dynamic_Active_Impl_i_q = iqPar;
                case 7
                    opGridAC.Loads_parallel_RL_static_Active_Impl_i_d = idPar;
                    opGridAC.Loads_parallel_RL_static_Active_Impl_i_q = iqPar;
                case 8
                    opGridAC.Loads_parallel_RL_static_omega_Active_Impl_i_d = idPar;
                    opGridAC.Loads_parallel_RL_static_omega_Active_Impl_i_q = iqPar;
            end
        end
    end
end
assignin(mdlWks, 'op_Grid_AC_Dynamic', opGridAC);

% -------------------------------------------------------------------------
% 3. DC network parameters, complexity vectors & operating point
% -------------------------------------------------------------------------
hasDC = ~isempty(TbusDC) && height(TbusDC) > 0;
DC = mdlWks.getVariable('DC');
opGridDC = mdlWks.getVariable('op_Grid_DC_Dynamic');

if hasDC
    DC.include = true;
    DC.Nbuses  = height(TbusDC);
    DC.bus     = TbusDC.bus_i(:);
    DC.v       = reshape(double(TbusDC.Vm), 1, []);
    gsDC       = double(TbusDC.Gs(:));
    DC.Rbus    = inf(size(gsDC));
    nzGs       = gsDC ~= 0;
    DC.Rbus(nzGs) = Sb ./ gsDC(nzGs);

    [G_dc, G_dc_noloads] = local_admitt(Sb, TbusDC, TlineDC);
    G_dc = real(G_dc);
    G_dc_noloads = real(G_dc_noloads);
    assignin(mdlWks, 'conductance_matrix', G_dc);
    assignin(mdlWks, 'iconductance_matrix', local_safe_inv(G_dc));
    assignin(mdlWks, 'conductance_matrix_noloads', G_dc_noloads);
    assignin(mdlWks, 'iconductance_matrix_noloads', local_safe_inv(G_dc_noloads));

    loadIdxDC = find(TbusDC.Pd ~= 0);
    if ~isempty(loadIdxDC)
        PloadsDC   = reshape(double(TbusDC.Pd(loadIdxDC)), 1, []);
        pSafe      = PloadsDC;
        pSafe(pSafe == 0) = 1e-10;
        R_Loads_DC = (DC.v(loadIdxDC).^2) .* Sb ./ pSafe;
        Load_busesDC = loadIdxDC(:);
        loadDC_complexity = local_encode_dc_load_complexity(TbusDC.Complexity_load(loadIdxDC));
    else
        Load_busesDC = 1;
        PloadsDC     = 0;
        R_Loads_DC   = 1e6;
        loadDC_complexity = 0;
    end
    assignin(mdlWks, 'Load_busesDC', Load_busesDC);
    assignin(mdlWks, 'PloadsDC', PloadsDC);
    assignin(mdlWks, 'R_Loads_DC', R_Loads_DC);
    assignin(mdlWks, 'loadDC_complexity', loadDC_complexity);

    hasDCLines = ~isempty(TlineDC) && height(TlineDC) > 0;
    if hasDCLines
        nLDC = height(TlineDC);
        dcFrom = zeros(nLDC, 1);
        dcTo   = zeros(nLDC, 1);
        for l = 1:nLDC
            dcFrom(l) = find(DC.bus == TlineDC.fbus(l), 1);
            dcTo(l)   = find(DC.bus == TlineDC.tbus(l), 1);
        end
        DC.from = dcFrom;
        DC.to   = dcTo;
        Mto   = zeros(DC.Nbuses, nLDC);
        Mfrom = zeros(DC.Nbuses, nLDC);
        for l = 1:nLDC
            Mto(dcTo(l), l)     = 1;
            Mfrom(dcFrom(l), l) = 1;
        end
        DC.Lines_toM   = sparse(Mto);
        DC.Lines_fromM = sparse(Mfrom);
        DC.R_Lines     = reshape(double(TlineDC.r), 1, []);
        DC.L_Lines     = reshape(double(TlineDC.x), 1, []);
        DC.CbusDC      = double(TbusDC.Bs(:)) / Sb;
        for b = 1:DC.Nbuses
            DC.CbusDC(b) = DC.CbusDC(b) + 0.5 * sum(TlineDC.b(dcFrom == b | dcTo == b));
        end
        lineDC_complexity = local_encode_dc_line_complexity(TlineDC.Complexity);
    else
        DC.from        = 1;
        DC.to          = 1;
        DC.Lines_toM   = sparse(DC.Nbuses, 1);
        DC.Lines_fromM = sparse(DC.Nbuses, 1);
        DC.R_Lines     = 1;
        DC.L_Lines     = 0.01;
        DC.CbusDC      = double(TbusDC.Bs(:)) / Sb;
        lineDC_complexity = 0;
    end

    % Apply .dyr DCBUS and DCLINE overrides before zero-inductance demotion
    if ~isempty(session.dyrData) && isfield(session.dyrData, 'records')
        for kr = 1:numel(session.dyrData.records)
            rc = session.dyrData.records(kr);
            if strcmp(rc.model, 'DCBUS')
                bIdx = find(DC.bus == rc.bus, 1);
                if ~isempty(bIdx)
                    if isfield(rc.params, 'C') && isnumeric(rc.params.C)
                        DC.CbusDC(bIdx) = double(rc.params.C(1));
                    end
                    if isfield(rc.params, 'R') && isnumeric(rc.params.R) && rc.params.R(1) > 0
                        DC.Rbus(bIdx) = double(rc.params.R(1));
                    end
                end
            elseif strcmp(rc.model, 'DCLINE') && isfield(rc.params, 'L') && isnumeric(rc.params.L) && hasDCLines
                lIdx = find((TlineDC.fbus == rc.fbus & TlineDC.tbus == rc.tbus) | ...
                            (TlineDC.fbus == rc.tbus & TlineDC.tbus == rc.fbus), 1);
                if ~isempty(lIdx), DC.L_Lines(lIdx) = double(rc.params.L(1)); end
            end
        end
    end

    if hasDCLines
        lineDC_complexity(lineDC_complexity == 1 & DC.L_Lines <= 0) = 2;
        DC.L_Lines(DC.L_Lines <= 0) = 0.01;
    end
    DC.CbusDC(DC.CbusDC <= 0) = 0.02;

    if hasDCLines
        DC.i_branch = (DC.v(DC.from) - DC.v(DC.to)) ./ max(DC.R_Lines, 1e-9);
    else
        DC.i_branch = 0;
    end

    busDC_complexity = local_encode_dc_bus_complexity(TbusDC.Complexity_bus);
    assignin(mdlWks, 'busDC_complexity', busDC_complexity);
    assignin(mdlWks, 'lineDC_complexity', lineDC_complexity);
    for kComp = 1:2
        assignin(mdlWks, sprintf('has_busDC_c%d', kComp), double(any(busDC_complexity == kComp)));
        assignin(mdlWks, sprintf('has_lineDC_c%d', kComp), double(hasDCLines && any(lineDC_complexity == kComp)));
    end
    for kComp = 1:4
        assignin(mdlWks, sprintf('has_loadDC_c%d', kComp), double(~isempty(loadIdxDC) && any(loadDC_complexity == kComp)));
    end

    allStaticDC = all(busDC_complexity == 2) && ...
                  all(loadDC_complexity == 2 | loadDC_complexity == 0) && ...
                  all(lineDC_complexity == 2 | lineDC_complexity == 0);
    if allStaticDC
        dc_grid_mode = 2;
    else
        dc_grid_mode = 1;
    end

    opGridDC.Buses_dynamic_Active_Impl_I1 = local_nonempty(DC.v(busDC_complexity == 1));
    opGridDC.Buses_static_Active_Impl_LPF_d1_I = local_nonempty(DC.v(busDC_complexity == 2));
    if hasDCLines && any(lineDC_complexity == 1)
        opGridDC.Lines_Dynamic_Active_Impl_Integrator = DC.i_branch(lineDC_complexity == 1);
    end
    opGridDC.Loads_DynamicP_Active_Impl_PI_Load_I = zeros(1, max(1, sum(loadDC_complexity == 4)));
else
    DC.include = false;
    dc_grid_mode = 0;
end
assignin(mdlWks, 'dc_grid_mode', dc_grid_mode);
assignin(mdlWks, 'op_Grid_DC_Dynamic', opGridDC);

% -------------------------------------------------------------------------
% 4. Consensus / Secondary & Tertiary agent arrays
% -------------------------------------------------------------------------
nACgen = height(TgenAC);
nDCgen = height(TgenDC);
nIC    = height(Tintconv);

buses_ac_gen = zeros(nACgen, 1);
Fb_ac_units  = ones(1, nACgen);
for g = 1:nACgen
    bPos = find(bus_ids == TgenAC.bus(g), 1);
    if isempty(bPos)
        error('VFlexP:BuildAndInit:BadGenBus', 'AC Generator %d references unknown AC bus %d.', g, TgenAC.bus(g));
    end
    buses_ac_gen(g) = bPos;
    Fb_ac_units(g)  = fref.buses(bPos);
end

Fb_ic_units = ones(1, nIC);
buses_ic_ac = ones(nIC, 1);
for k = 1:nIC
    bPos = find(bus_ids == Tintconv.tbus(k), 1);
    if isempty(bPos)
        error('VFlexP:BuildAndInit:BadICBus', 'Converter %d references unknown AC tbus %d.', k, Tintconv.tbus(k));
    end
    buses_ic_ac(k) = bPos;
    Fb_ic_units(k) = fref.buses(bPos);
end
if nDCgen == 0
    Fb_dc_units = [];
elseif nIC > 0
    Fb_dc_units = repmat(Fb_ic_units(1), 1, nDCgen);
else
    Fb_dc_units = ones(1, nDCgen);
end
fref.Fb_units = [Fb_ac_units, Fb_dc_units, Fb_ic_units];

Ps_ac = reshape(double(TgenAC.Pg), 1, []) / Sb;
Qs_ac = reshape(double(TgenAC.Qg), 1, []) / Sb;
vs_ac = v_buses_mod(buses_ac_gen);

buses_lf = buses_ac_gen;
Ps_lf    = Ps_ac;
Qs_lf    = Qs_ac;
vs_lf    = vs_ac;
is_dc_agent      = false(1, nACgen);
is_intconv_agent = false(1, nACgen);

if nDCgen > 0
    Ps_dc = reshape(double(TgenDC.Pg), 1, []) / Sb;
    vdc_ref = ones(1, nDCgen);
    for g = 1:nDCgen
        bPos = find(DC.bus == TgenDC.bus(g), 1);
        if ~isempty(bPos), vdc_ref(g) = DC.v(bPos); end
    end
    if nIC > 0
        acLinkBus = buses_ic_ac(1);
    elseif nACgen > 0
        acLinkBus = buses_ac_gen(1);
    else
        acLinkBus = 1;
    end
    buses_lf = [buses_lf; repmat(acLinkBus, nDCgen, 1)];
    Ps_lf    = [Ps_lf, Ps_dc];
    Qs_lf    = [Qs_lf, zeros(1, nDCgen)];
    vs_lf    = [vs_lf, vdc_ref];
    is_dc_agent      = [is_dc_agent, true(1, nDCgen)];
    is_intconv_agent = [is_intconv_agent, false(1, nDCgen)];
end

intconv_PF = session.pfResult.intconv_PF;
assignin(mdlWks, 'intconv_PF', intconv_PF);

if nIC > 0
    if isfield(intconv_PF, 'include') && intconv_PF.include
        Ps_ic = -reshape(double(intconv_PF.P), 1, []);
        Qs_ic = -reshape(double(intconv_PF.Q), 1, []);
        vs_ic =  reshape(double(intconv_PF.vmod), 1, []);
    else
        Ps_ic = zeros(1, nIC);
        Qs_ic = zeros(1, nIC);
        vs_ic = v_buses_mod(buses_ic_ac);
    end
    buses_lf = [buses_lf; buses_ic_ac];
    Ps_lf    = [Ps_lf, Ps_ic];
    Qs_lf    = [Qs_lf, Qs_ic];
    vs_lf    = [vs_lf, vs_ic];
    is_dc_agent      = [is_dc_agent, false(1, nIC)];
    is_intconv_agent = [is_intconv_agent, true(1, nIC)];
end

if isempty(buses_lf)
    buses_lf = 1;
    Ps_lf = 1e-18;
    Qs_lf = 1e-18;
    vs_lf = 1;
    is_dc_agent = false;
    is_intconv_agent = false;
    fref.Fb_units = 1;
end

Ps_safe = Ps_lf;
Ps_safe(Ps_safe == 0) = 1e-18;
Qs_safe = Qs_lf;
Qs_safe(~is_dc_agent & Qs_safe == 0) = 1e-18;
iPs_lf = 1 ./ Ps_safe;
iPs_lf(~isfinite(iPs_lf)) = 1e18;
iQs_lf = zeros(size(Qs_safe));
iQs_lf(~is_dc_agent) = 1 ./ Qs_safe(~is_dc_agent);
iQs_lf(~isfinite(iQs_lf)) = 1e18;

MAS_Pnom = reshape(double(TgenAC.mBase), 1, []) / Sb;
if nDCgen > 0
    MAS_Pnom = [MAS_Pnom, reshape(double(TgenDC.mBase), 1, []) / Sb];
end
if nIC > 0
    sNomIC = reshape(double(Tintconv.rateA), 1, []);
    sNomIC(sNomIC <= 0) = Sb;
    MAS_Pnom = [MAS_Pnom, sNomIC / Sb];
end
if isempty(MAS_Pnom), MAS_Pnom = 1; end

gen_static_params = struct('mBase', double(TgenAC.mBase(:)));
if isempty(gen_static_params.mBase), gen_static_params.mBase = Sb; end

assignin(mdlWks, 'gen_static_params', gen_static_params);
assignin(mdlWks, 'DGs', cellstr(TgenAC.type));
assignin(mdlWks, 'PQ_gen', TbusAC.type(buses_ac_gen));
assignin(mdlWks, 'buses_lf', buses_lf);
assignin(mdlWks, 'Ps_lf', Ps_lf);
assignin(mdlWks, 'Qs_lf', Qs_lf);
assignin(mdlWks, 'iPs_lf', iPs_lf);
assignin(mdlWks, 'iQs_lf', iQs_lf);
assignin(mdlWks, 'vs_lf', vs_lf);
assignin(mdlWks, 'is_dc_agent', is_dc_agent);
assignin(mdlWks, 'is_intconv_agent', is_intconv_agent);
assignin(mdlWks, 'MAS_Pnom', MAS_Pnom);
assignin(mdlWks, 'sc_mode', local_encode_sc_mode(session.config.sim.sc_structure));
assignin(mdlWks, 'tc_mode', local_encode_tc_mode(session.config.sim.tc_structure));

% -------------------------------------------------------------------------
% 5. Configure parameters & initialize AC Generators, DC Generators, IntConvs
% -------------------------------------------------------------------------
warningsList = strings(0, 1);
unitCounts   = struct();
opUnits      = struct();

acUnitNames = {'DFIG','GFL','GFM','IG','MGES','MSG','SC','SG','VSM'};
for k = 1:numel(acUnitNames)
    uName = acUnitNames{k};
    rowsU = find(strcmp(TgenAC.type, uName));
    nU    = numel(rowsU);
    unitCounts.(uName) = nU;
    if nU == 0
        assignin(mdlWks, ['enable_' uName], 0);
        fref.(uName) = 1;
        continue;
    end

    assignin(mdlWks, ['enable_' uName], 1);
    uBuses = buses_ac_gen(rowsU);
    uToM   = zeros(Nbuses, nU);
    uToM(sub2ind([Nbuses, nU], uBuses(:)', 1:nU)) = 1;
    fref.(uName) = fref.buses(uBuses(:)');

    assignin(mdlWks, [uName 's_buses'], uBuses(:));
    assignin(mdlWks, [uName '_idx'], rowsU(:));
    assignin(mdlWks, [uName 's_toM'], uToM);

    prmU = local_build_ac_unit_prm(uName, rowsU, TgenAC, TbusAC, TFb, Sb, session, manifest);
    assignin(mdlWks, ['prm_' uName], prmU);
end

% Assign fref before unit initialization harnesses read it
dcUnitNames = {'DCBUCK', 'DC_PES'};
for k = 1:numel(dcUnitNames)
    uName = dcUnitNames{k};
    fref.(uName) = 1;
end
icUnitNames = {'intconv_GFL', 'intconv_GFM', 'intconv_VSM'};
for k = 1:numel(icUnitNames)
    uName = icUnitNames{k};
    rowsIC = find(strcmp(Tintconv.type, uName));
    if ~isempty(rowsIC)
        fref.(uName) = fref.buses(buses_ic_ac(rowsIC(:))');
    else
        fref.(uName) = 1;
    end
end
assignin(mdlWks, 'fref', fref);

% Initialize active AC generators
for k = 1:numel(acUnitNames)
    uName = acUnitNames{k};
    rowsU = find(strcmp(TgenAC.type, uName));
    if isempty(rowsU), continue; end

    uBuses = buses_ac_gen(rowsU);
    vU     = v_buses_mod(uBuses(:)');
    thU    = v_buses_angle(uBuses(:)');
    pU     = Ps_ac(rowsU(:)');
    qU     = Qs_ac(rowsU(:)');
    prmU   = mdlWks.getVariable(['prm_' uName]);
    opU0   = mdlWks.getVariable(['op_' uName]);

    if strcmp(uName, 'GFM') && ismember(lower(strtrim(char(string(prmU.dynamic_VCC)))), {'static', 'dynamic'})
        opU = local_init_op_GFM(opU0, prmU, vU, thU, pU, qU);
    elseif strcmp(uName, 'GFL')
        opU = local_init_op_GFL(opU0, prmU, vU, thU, pU, qU);
    else
        [opU, warnMsg] = local_findop_ac_unit(uName, manifest.units.(uName), prmU, opU0, ...
            vU, thU, pU, qU, mdlWks, omega_base, Sb);
        if strlength(warnMsg) > 0
            warningsList(end+1, 1) = warnMsg; %#ok<AGROW>
        end
    end
    assignin(mdlWks, ['op_' uName], opU);
    opUnits.(uName) = opU;
end

% Configure & initialize DC Generators
for k = 1:numel(dcUnitNames)
    uName = dcUnitNames{k};
    rowsD = find(strcmp(TgenDC.type, uName));
    nD    = numel(rowsD);
    unitCounts.(uName) = nD;
    if nD == 0 || ~hasDC
        assignin(mdlWks, ['enable_' uName], 0);
        continue;
    end

    assignin(mdlWks, ['enable_' uName], 1);
    dcBusesU = zeros(1, nD);
    for ii = 1:nD
        dcBusesU(ii) = find(DC.bus == TgenDC.bus(rowsD(ii)), 1);
    end
    dcToMU = zeros(DC.Nbuses, nD);
    dcToMU(sub2ind([DC.Nbuses, nD], dcBusesU, 1:nD)) = 1;

    DC.([uName 's_buses'])     = dcBusesU;
    DC.([uName 's_toM'])       = dcToMU;
    DC.([uName 's_Pg'])        = reshape(double(TgenDC.Pg(rowsD)), 1, []);
    DC.([uName 's_agent_idx']) = nACgen + reshape(rowsD, 1, []);
    assignin(mdlWks, 'DC', DC);

    prmD = local_build_dc_unit_prm(uName, rowsD, TgenDC, Sb, session, manifest);
    assignin(mdlWks, ['prm_' uName], prmD);

    opD0 = mdlWks.getVariable(['op_' uName]);
    [opD, warnMsg] = local_findop_dc_unit(uName, manifest.units.(uName), prmD, opD0, DC, mdlWks, Sb);
    if strlength(warnMsg) > 0
        warningsList(end+1, 1) = warnMsg; %#ok<AGROW>
    end
    assignin(mdlWks, ['op_' uName], opD);
    opUnits.(uName) = opD;
end

% Configure & initialize DC/AC Interlinking Converters
for k = 1:numel(icUnitNames)
    uName  = icUnitNames{k};
    rowsIC = find(strcmp(Tintconv.type, uName));
    nC     = numel(rowsIC);
    unitCounts.(uName) = nC;
    if nC == 0 || ~hasDC
        assignin(mdlWks, ['enable_' uName], 0);
        continue;
    end

    assignin(mdlWks, ['enable_' uName], 1);
    acBusesC = buses_ic_ac(rowsIC);
    acToMC   = zeros(Nbuses, nC);
    acToMC(sub2ind([Nbuses, nC], acBusesC(:)', 1:nC)) = 1;

    dcBusesC = zeros(nC, 1);
    for ii = 1:nC
        dcBusesC(ii) = find(DC.bus == Tintconv.fbus(rowsIC(ii)), 1);
    end
    dcToMC = zeros(DC.Nbuses, nC);
    dcToMC(sub2ind([DC.Nbuses, nC], dcBusesC(:)', 1:nC)) = 1;

    DC.([uName 's_buses'])     = acBusesC(:);
    DC.([uName 's_DC_buses'])  = dcBusesC(:);
    DC.([uName 's_toM'])       = acToMC;
    DC.([uName 's_DC_toM'])    = dcToMC;
    DC.([uName 's_agent_idx']) = nACgen + nDCgen + reshape(rowsIC, 1, []);

    dparams = struct();
    dparams.ratio  = reshape(double(Tintconv.ratio(rowsIC)), 1, []);
    dparams.ratio(dparams.ratio == 0) = 1.0;
    dparams.type   = reshape(double(Tintconv.CONV(rowsIC)), 1, []);
    dparams.Rf     = reshape(double(Tintconv.r(rowsIC)), 1, []);
    dparams.Lf     = reshape(double(Tintconv.x(rowsIC)), 1, []);
    dparams.alpha1 = reshape(double(Tintconv.ALPHA1(rowsIC)), 1, []);
    dparams.alpha2 = reshape(double(Tintconv.ALPHA2(rowsIC)), 1, []);
    dparams.alpha3 = reshape(double(Tintconv.ALPHA3(rowsIC)), 1, []);
    dparams.Gsw    = reshape(double(Tintconv.GSW(rowsIC)), 1, []);

    prmC = local_build_ic_unit_prm(uName, rowsIC, Tintconv, TFb, Sb, session, manifest);
    assignin(mdlWks, ['prm_' uName], prmC);

    if ismember(uName, {'intconv_VSM', 'intconv_GFM'})
        if isfield(prmC, 'R_out') && numel(prmC.R_out) == nC
            dparams.Rf = reshape(double(prmC.R_out), 1, []);
        end
        if isfield(prmC, 'L_out') && numel(prmC.L_out) == nC
            dparams.Lf = reshape(double(prmC.L_out), 1, []);
        end
    end

    vC  = reshape(double(intconv_PF.vmod(rowsIC)), 1, []);
    thC = deg2rad(reshape(double(intconv_PF.vangle(rowsIC)), 1, []));
    pC  = -reshape(double(intconv_PF.P(rowsIC)), 1, []);
    qC  = -reshape(double(intconv_PF.Q(rowsIC)), 1, []);

    if isfield(intconv_PF, 'Pdc') && ~isempty(intconv_PF.Pdc)
        if strcmp(uName, 'intconv_GFM')
            hasPoly = (dparams.alpha1 ~= 0) | (dparams.alpha2 ~= 0) | (dparams.alpha3 ~= 0);
            if any(hasPoly)
                iMag    = hypot(pC, qC) ./ max(vC, 1e-6);
                pBridge = pC + dparams.Rf .* (iMag.^2);
                pDcReq  = reshape(double(intconv_PF.Pdc(rowsIC)), 1, []);
                pDcMod  = pBridge + dparams.alpha1 + dparams.alpha2 .* iMag + dparams.alpha3 .* (iMag.^2);
                dparams.alpha1(hasPoly) = dparams.alpha1(hasPoly) + (pDcReq(hasPoly) - pDcMod(hasPoly));
            end
        else
            pDcReq = reshape(double(intconv_PF.Pdc(rowsIC)), 1, []);
            pDcMod = pC + dparams.alpha1 + dparams.alpha2 .* hypot(pC, 1e-6) + dparams.alpha3 .* (pC.^2);
            dparams.alpha1 = dparams.alpha1 + (pDcReq - pDcMod);
            dparams.Gsw(:) = 0;
        end
    end

    DC.([uName '_dparams']) = dparams;
    assignin(mdlWks, 'DC', DC);

    opC0 = mdlWks.getVariable(['op_' uName]);
    if strcmp(uName, 'intconv_GFL')
        opC = local_init_op_intconv_GFL(opC0, prmC, dparams, vC, thC, pC, qC);
    elseif strcmp(uName, 'intconv_GFM')
        opC = local_init_op_intconv_GFM(opC0, prmC, vC, thC, pC, qC);
    else
        opC = local_init_op_intconv_VSM(opC0, prmC, dparams, vC, thC, pC, qC);
    end
    assignin(mdlWks, ['op_' uName], opC);
    opUnits.(uName) = opC;
end

assignin(mdlWks, 'DC', DC);

% Update diagram to propagate active variants & verify clean compilation
set_param('VFlexP', 'SimulationCommand', 'update');
x0Struct = Simulink.BlockDiagram.getInitialState('VFlexP');
nStates = 0;
if isstruct(x0Struct) && isfield(x0Struct, 'signals')
    for s = 1:numel(x0Struct.signals)
        nStates = nStates + numel(x0Struct.signals(s).values);
    end
end

elapsed = toc(tStart);
initRes = struct( ...
    'ok', isempty(warningsList), ...
    'nStates', nStates, ...
    'unitCounts', unitCounts, ...
    'acGridMode', ac_grid_mode, ...
    'dcGridMode', dc_grid_mode, ...
    'opUnits', opUnits, ...
    'opGridAC', opGridAC, ...
    'opGridDC', opGridDC, ...
    'warnings', warningsList, ...
    'elapsed', elapsed);

session.initResult = initRes;
session.log('Initialized dynamic model in %.3f s (%d active states).', elapsed, nStates);
for w = 1:numel(warningsList)
    session.log('WARNING: %s', warningsList(w));
end
end

% =========================================================================
% Local Helpers
% =========================================================================
function [ref_units, use_sync_frame] = local_resolve_area_frames(areas, TbusAC, TgenAC, session)
nAreas = numel(areas);
ref_units      = ones(1, nAreas);
use_sync_frame = ones(1, nAreas);

policy = 'first-gen';
if isfield(session.config, 'model') && isfield(session.config.model, 'frame_policy') && ...
        ~isempty(session.config.model.frame_policy)
    policy = lower(strtrim(char(string(session.config.model.frame_policy))));
end

if ismember(policy, {'fixed', 'nominal', 'sync', 'coi'})
    return;
end

if isempty(TgenAC) || height(TgenAC) == 0
    return;
end

gfmTypes = {'MGES', 'SG', 'MSG', 'SC', 'GFM', 'VSM'};
isGFMGen = ismember(cellstr(string(TgenAC.type)), gfmTypes);
if ~isempty(session.dyrData) && isstruct(session.dyrData) && isfield(session.dyrData, 'records')
    recs = session.dyrData.records;
    for g = 1:height(TgenAC)
        if strcmp(TgenAC.type(g), "VSM")
            r = VFlexPSession.findBusRecord(recs, round(TgenAC.bus(g)), {'device'});
            if ~isempty(r) && isfield(r.params, 'Hv') && isnumeric(r.params.Hv) && ~isfinite(r.params.Hv(1))
                isGFMGen(g) = false;
            end
        end
    end
end
isGFMGen = isGFMGen(:);

unitTargetBus = NaN;
if startsWith(policy, 'unit:')
    tok = strtrim(extractAfter(policy, 'unit:'));
    unitTargetBus = str2double(tok);
end

for a = 1:nAreas
    inArea = TbusAC.bus_i(TbusAC.area == areas(a));
    if isempty(inArea)
        continue;
    end
    gPos = [];
    if ~isnan(unitTargetBus) && ismember(unitTargetBus, inArea)
        gPos = find(TgenAC.bus == unitTargetBus & isGFMGen, 1);
        if isempty(gPos)
            gPos = find(TgenAC.bus == unitTargetBus, 1);
        end
    elseif strcmp(policy, 'slack')
        slackBuses = TbusAC.bus_i(TbusAC.area == areas(a) & TbusAC.type == 3);
        gPos = find(ismember(TgenAC.bus, slackBuses) & isGFMGen, 1);
        if isempty(gPos)
            gPos = find(ismember(TgenAC.bus, inArea) & isGFMGen, 1);
        end
    else
        % 'first-gen' (or non-target area in multi-area 'unit:<id>')
        gPos = find(ismember(TgenAC.bus, inArea) & isGFMGen, 1);
    end
    if ~isempty(gPos)
        ref_units(a)      = gPos;
        use_sync_frame(a) = 0;
    end
end
end

function manifest = local_get_manifest(repoRoot)
persistent cachedManifest cachedPath
mPath = fullfile(repoRoot, 'Simulink', 'lib', 'vflexp_lib_manifest.mat');
if isempty(cachedManifest) || ~strcmp(cachedPath, mPath)
    S = load(mPath, 'manifest');
    cachedManifest = S.manifest;
    cachedPath = mPath;
end
manifest = cachedManifest;
end

function out = local_nonempty(v)
if isempty(v)
    out = 0;
else
    out = reshape(double(v), 1, []);
end
end

function invM = local_safe_inv(M)
if isempty(M)
    invM = 1;
elseif rcond(full(M)) < 1e-13
    invM = pinv(full(M));
else
    invM = inv(full(M));
end
end

function [Y, Ynoloads] = local_admitt(Sb, Tbus, Tline)
nbus = height(Tbus);
busn = Tbus.bus_i;
aux = zeros(nbus, 1);
if ~isempty(Tline) && height(Tline) > 0
    for i = 1:nbus
        aux(i) = 0.5 * sum(Tline.b(Tline.fbus == busn(i) | Tline.tbus == busn(i)));
    end
end
gs = double(Tbus.Gs);
Cbus = double(Tbus.Bs) / Sb + aux;
hasLd = (Tbus.Pd ~= 0) | (Tbus.Qd ~= 0);
Yloads = zeros(nbus, 1);
if any(hasLd)
    Pl = double(Tbus.Pd(hasLd)) / Sb;
    Ql = double(Tbus.Qd(hasLd)) / Sb;
    vm = double(Tbus.Vm(hasLd));
    s2 = max(Pl.^2 + Ql.^2, 1e-20);
    RL = (vm.^2) .* (Pl ./ s2);
    LL = (vm.^2) .* (Ql ./ s2);
    Yloads(hasLd) = 1 ./ (RL + 1i * LL);
end

Y = zeros(nbus, nbus);
if ~isempty(Tline) && height(Tline) > 0
    tap = double(Tline.ratio);
    tap(tap == 0) = 1.0;
    ph = deg2rad(double(Tline.angle));
    r = double(Tline.r);
    x = double(Tline.x);
    for j = 1:height(Tline)
        i1 = find(busn == Tline.fbus(j), 1);
        i2 = find(busn == Tline.tbus(j), 1);
        yij = 1 / (r(j) + 1i * x(j));
        Y(i1, i1) = Y(i1, i1) + yij / (tap(j)^2);
        Y(i2, i2) = Y(i2, i2) + yij;
        Y(i1, i2) = Y(i1, i2) - yij / (tap(j) * (cos(ph(j)) - 1i * sin(ph(j))));
        Y(i2, i1) = Y(i2, i1) - yij / (tap(j) * (cos(ph(j)) + 1i * sin(ph(j))));
    end
end
for j = 1:nbus
    Y(j, j) = Y(j, j) + Yloads(j) + (gs(j) / Sb) + 1i * Cbus(j);
end
Ynoloads = Y - diag(Yloads);
end

function c = local_encode_ac_bus_complexity(strs)
n = numel(strs);
c = 3 * ones(1, n);
opts = {'dynamic', 'static', 'static omega', 'only r'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if strcmp(s, 'emt'), s = 'dynamic';
    elseif strcmp(s, 'rms'), s = 'static omega';
    end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_encode_ac_line_complexity(strs)
n = numel(strs);
c = 3 * ones(1, n);
opts = {'dynamic', 'static', 'static omega'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if strcmp(s, 'emt'), s = 'dynamic';
    elseif strcmp(s, 'rms'), s = 'static omega';
    end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_encode_ac_load_complexity(strs)
n = numel(strs);
c = 3 * ones(1, n);
opts = {'series rl dynamic', 'series rl static', 'series rl static omega', ...
        'pq', 'dynamic pq', 'parallel rl dynamic', 'parallel rl static', ...
        'parallel rl static omega', 'constant current'};
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if strcmp(s, 'emt') || strcmp(s, 'dynamic')
        s = 'series rl dynamic';
    elseif strcmp(s, 'rms') || strcmp(s, 'static omega')
        s = 'series rl static omega';
    elseif strcmp(s, 'static')
        s = 'series rl static';
    end
    idx = find(strcmp(s, opts), 1);
    if ~isempty(idx), c(k) = idx; end
end
end

function c = local_encode_dc_bus_complexity(strs)
n = numel(strs);
c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc dynamic', 'dc emt', 'dynamic', 'emt'})
        c(k) = 1;
    end
end
end

function c = local_encode_dc_line_complexity(strs)
n = numel(strs);
c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc dynamic', 'dc emt', 'dynamic', 'emt'})
        c(k) = 1;
    end
end
end

function c = local_encode_dc_load_complexity(strs)
n = numel(strs);
c = 2 * ones(1, n);
for k = 1:n
    s = lower(strtrim(char(string(strs(k)))));
    if ismember(s, {'dc constr', 'constr'})
        c(k) = 1;
    elseif ismember(s, {'dc constp', 'constp', 'dc static', 'dc static omega', 'dc rms'})
        c(k) = 2;
    elseif ismember(s, {'dc consti', 'consti'})
        c(k) = 3;
    elseif ismember(s, {'dc dynamicp', 'dynamicp', 'dc dynamic', 'dc emt'})
        c(k) = 4;
    end
end
end

function m = local_encode_sc_mode(scName)
s = lower(strtrim(char(string(scName))));
switch s
    case 'pilot nodes', m = 1;
    case 'matrices', m = 2;
    case 'matrices1', m = 3;
    case 'matrices with delay', m = 4;
    case 'matrices with delay and limitation', m = 5;
    case 'mas messages', m = 6;
    otherwise, m = 0;
end
end

function m = local_encode_tc_mode(tcName)
s = lower(strtrim(char(string(tcName))));
switch s
    case 'scheduled', m = 1;
    case 'mas', m = 2;
    case 'mas edp', m = 3;
    case 'mas messages', m = 4;
    otherwise, m = 0;
end
end

% -------------------------------------------------------------------------
% Unit Parameter Builders (.dyr + MBASE scaling + unitOverrides)
% -------------------------------------------------------------------------
function prm = local_build_ac_unit_prm(uName, rowsU, TgenAC, TbusAC, TFb, Sb, session, manifest)
prm = manifest.units.(uName).defaultPrm;
nU  = numel(rowsU);
hasDyr = ~isempty(session.dyrData) && isstruct(session.dyrData) && isfield(session.dyrData, 'records');
if hasDyr
    recs = session.dyrData.records;
else
    recs = [];
end

switch uName
    case {'MGES', 'SG'}
        prm.SGs_nodes = 1:nU;
        if hasDyr
            machRecs = cell(1, nU);
            hasMach  = false;
            for u = 1:nU
                bId = round(TgenAC.bus(rowsU(u)));
                machRecs{u} = VFlexPSession.findBusRecord(recs, bId, {'machine'});
                if ~isempty(machRecs{u}), hasMach = true; end
            end
            if hasMach
                Tpdo  = 8.0 * ones(1, nU);  Tppdo = 0.03 * ones(1, nU);
                Tpqo  = 0.4 * ones(1, nU);  Tppqo = 0.05 * ones(1, nU);
                H     = 6.5 * ones(1, nU);  D     = zeros(1, nU);
                Xd    = 1.8 * ones(1, nU);  Xq    = 1.7 * ones(1, nU);
                Xpd   = 0.3 * ones(1, nU);  Xpq   = 0.55 * ones(1, nU);
                Xppd  = 0.25 * ones(1, nU); Xppq  = 0.25 * ones(1, nU);
                Xl    = 0.2 * ones(1, nU);  Ra    = 0.0025 * ones(1, nU);
                fBase = double(TFb.Fb(1)) * ones(1, nU);
                convMode = 'classical';
                detailMode = 'EMT8';

                ovNames = {'Lfd','Rfd','Lkd1','Rkd1','Lkq1','Rkq1','Lkq2','Rkq2'};
                ovStruct = struct();
                for oi = 1:numel(ovNames)
                    onm = ovNames{oi};
                    vals = nan(1, nU);
                    for u = 1:nU
                        if ~isempty(machRecs{u})
                            vals(u) = local_pget(machRecs{u}.params, onm, NaN);
                        end
                    end
                    if all(~isnan(vals)), ovStruct.(onm) = vals; end
                end

                for u = 1:nU
                    bId = round(TgenAC.bus(rowsU(u)));
                    bIdx = find(TbusAC.bus_i == bId, 1);
                    if ~isempty(bIdx)
                        fIdx = find(TFb.Area == TbusAC.area(bIdx), 1);
                        if ~isempty(fIdx), fBase(u) = double(TFb.Fb(fIdx)); end
                    end
                    mr = machRecs{u};
                    if ~isempty(mr)
                        p = mr.params;
                        Tpdo(u)  = local_pget(p, 'Tpdo', 8.0);
                        Tppdo(u) = local_pget(p, 'Tppdo', 0.03);
                        Tppqo(u) = local_pget(p, 'Tppqo', 0.05);
                        H(u)     = local_pget(p, 'H', 6.5);
                        D(u)     = local_pget(p, 'D', 0.0);
                        Xd(u)    = local_pget(p, 'Xd', 1.8);
                        Xq(u)    = local_pget(p, 'Xq', 1.7);
                        Xpd(u)   = local_pget(p, 'Xpd', 0.3);
                        Xppd(u)  = local_pget(p, 'Xppd', 0.25);
                        Xppq(u)  = local_pget(p, 'XPPQ', Xppd(u));
                        Xl(u)    = local_pget(p, 'Xl', 0.2);
                        Ra(u)    = local_pget(p, 'RA', local_pget(p, 'Ra', 0.0));
                        if strcmp(mr.model, 'GENSAL')
                            Tpqo(u) = local_pget(p, 'Tpqo', max(Tppqo(u) * 10, 0.4));
                            Xpq(u)  = local_pget(p, 'Xpq', 0.5 * (Xq(u) + Xppq(u)));
                        else
                            Tpqo(u) = local_pget(p, 'Tpqo', 0.4);
                            Xpq(u)  = local_pget(p, 'Xpq', 0.55);
                        end
                        if isfield(p, 'CONVERSION')
                            convMode = lower(char(string(p.CONVERSION)));
                        end
                        if isfield(mr.complexity, 'MODE')
                            detailMode = upper(char(string(mr.complexity.MODE)));
                        elseif isfield(mr.complexity, 'MODEL')
                            detailMode = upper(char(string(mr.complexity.MODEL)));
                        end
                    end
                end

                fnDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions');
                oldPwd = pwd;
                pwdClean = onCleanup(@() cd(oldPwd));
                cd(fnDir);
                circ = param_genrou(Tpdo, Tppdo, Tpqo, Tppqo, H, D, Xd, Xq, Xpd, Xpq, Xppd, Xl, ...
                    Xppq, Ra, 0.0, fBase, convMode, ovStruct);
                clear pwdClean;
                prm.FBASE   = fBase(1);
                for fn = {'Rs_pu','Lmd_pu','Lmq_pu','Lssd_pu','Lssq_pu','Lfd_pu','Rfd_pu', ...
                          'Lkd1_pu','Rkd1_pu','Lkd2_pu','Rkd2_pu','Lkq1_pu','Rkq1_pu', ...
                          'Lkq2_pu','Rkq2_pu','Lkq3_pu','Rkq3_pu','H','D'}
                    prm.(fn{1}) = circ.(fn{1});
                end
                if strcmp(detailMode, 'RMS4')
                    prm.dynamic_stator = 'order 4';
                elseif strcmp(detailMode, 'CLASSICAL2')
                    prm.dynamic_stator = 'order 2';
                else
                    prm.dynamic_stator = 'order 8';
                end

                if strcmp(uName, 'MGES')
                    % AVR (ST4B / ST1)
                    prm.T_v  = 0.02 * ones(1, nU); prm.K_PR = zeros(1, nU); prm.K_IR = zeros(1, nU);
                    prm.K_A  = ones(1, nU);        prm.T_e  = 0.02 * ones(1, nU);
                    prm.T_1  = zeros(1, nU);       prm.T_2  = zeros(1, nU);
                    prm.K_P  = zeros(1, nU);       prm.K_C  = zeros(1, nU);
                    prm.K_PM = 1;                  prm.K_IM = 0;
                    for u = 1:nU
                        bId = round(TgenAC.bus(rowsU(u)));
                        ar = VFlexPSession.findBusRecord(recs, bId, {'avr'});
                        if ~isempty(ar) && ~strcmpi(local_cget(ar.complexity, 'MODE', 'ON'), 'OFF')
                            p = ar.params;
                            prm.T_v(u)  = max(local_pget(p, 'T_v', local_pget(p, 'T_r', 0.02)), 1e-3);
                            prm.K_PR(u) = local_pget(p, 'K_PR', 10.0);
                            prm.K_IR(u) = local_pget(p, 'K_IR', 10.0);
                            prm.K_A(u)  = local_pget(p, 'K_A', local_pget(p, 'K_PM', 1.0));
                            prm.T_e(u)  = max(local_pget(p, 'T_e', local_pget(p, 'T_A', 0.02)), 1e-3);
                            prm.T_1(u)  = local_pget(p, 'T_1', local_pget(p, 'T_B', 0.0));
                            prm.T_2(u)  = local_pget(p, 'T_2', local_pget(p, 'T_C', 0.0));
                            prm.K_P(u)  = local_pget(p, 'K_P', 0.0);
                            prm.K_C(u)  = local_pget(p, 'K_C', 0.0);
                        end
                    end

                    % Governor (IEEEG1)
                    prm.R_p = 0.05 * ones(1, nU); prm.T_3 = 0.1 * ones(1, nU);
                    prm.T_4 = 0.3 * ones(1, nU);  prm.T_5 = 5.0 * ones(1, nU);
                    prm.T_6 = 0.5 * ones(1, nU);  prm.K_1 = zeros(1, nU);
                    prm.K_3 = zeros(1, nU);       prm.K_5 = zeros(1, nU);
                    prm.T_meas = 0.01 * ones(1, nU); prm.cosn = ones(1, nU);
                    for u = 1:nU
                        bId = round(TgenAC.bus(rowsU(u)));
                        gr = VFlexPSession.findBusRecord(recs, bId, {'gov'});
                        if ~isempty(gr) && ~strcmpi(local_cget(gr.complexity, 'MODE', 'ON'), 'OFF')
                            p = gr.params;
                            prm.R_p(u)    = max(local_pget(p, 'R_p', local_pget(p, 'K', 0.05)), 1e-6);
                            prm.T_3(u)    = max(local_pget(p, 'T_3', 0.1), 1e-3);
                            prm.T_4(u)    = max(local_pget(p, 'T_4', 0.3), 1e-3);
                            prm.T_5(u)    = max(local_pget(p, 'T_5', 5.0), 1e-3);
                            prm.T_6(u)    = max(local_pget(p, 'T_6', 0.5), 1e-3);
                            prm.K_1(u)    = local_pget(p, 'K_1', 0.3);
                            prm.K_3(u)    = local_pget(p, 'K_3', 0.3);
                            prm.K_5(u)    = local_pget(p, 'K_5', 0.4);
                            prm.T_meas(u) = local_pget(p, 'T_meas', 0.01);
                            prm.cosn(u)   = local_pget(p, 'cosn', 1.0);
                        end
                    end

                    % Stabiliser (PSS2B)
                    prm.T_w1 = 2.0 * ones(1, nU);  prm.T_w2 = 2.0 * ones(1, nU); prm.T_w3 = 2.0 * ones(1, nU);
                    prm.T_7p = 2.0 * ones(1, nU);  prm.T_9  = 0.1 * ones(1, nU); prm.M = 5;
                    prm.K_S1 = zeros(1, nU);       prm.K_S2 = 0.15 * ones(1, nU); prm.K_S3 = ones(1, nU);
                    prm.T_1p = 0.15 * ones(1, nU); prm.T_2p = 0.03 * ones(1, nU);
                    prm.T_3p = 0.15 * ones(1, nU); prm.T_4p = 0.03 * ones(1, nU);
                    for u = 1:nU
                        bId = round(TgenAC.bus(rowsU(u)));
                        pr = VFlexPSession.findBusRecord(recs, bId, {'pss'});
                        if ~isempty(pr) && ~strcmpi(local_cget(pr.complexity, 'MODE', 'ON'), 'OFF')
                            p = pr.params;
                            ks1 = local_pget(p, 'K_S1', 0.0);
                            if ks1 ~= 0.0
                                prm.T_w1(u) = max(local_pget(p, 'T_w1', 2.0), 1e-3);
                                prm.T_w2(u) = max(local_pget(p, 'T_w2', 2.0), 1e-3);
                                prm.T_w3(u) = max(local_pget(p, 'T_w3', prm.T_w1(u)), 1e-3);
                                prm.T_7p(u) = max(local_pget(p, 'T_7p', local_pget(p, 'T_7', 2.0)), 1e-3);
                                prm.T_9(u)  = max(local_pget(p, 'T_9', 0.1), 1e-3);
                                prm.M       = max(round(local_pget(p, 'M', 5)), 0);
                                prm.K_S1(u) = ks1;
                                prm.K_S2(u) = local_pget(p, 'K_S2', 0.15);
                                prm.K_S3(u) = local_pget(p, 'K_S3', 1.0);
                                prm.T_1p(u) = local_pget(p, 'T_1p', local_pget(p, 'T_1', 0.15));
                                prm.T_2p(u) = max(local_pget(p, 'T_2p', local_pget(p, 'T_2', 0.03)), 1e-3);
                                prm.T_3p(u) = local_pget(p, 'T_3p', local_pget(p, 'T_3', 0.15));
                                prm.T_4p(u) = max(local_pget(p, 'T_4p', local_pget(p, 'T_4', 0.03)), 1e-3);
                            end
                        end
                    end
                end
            end
        end

    case 'VSM'
        if hasDyr
            vsmRecs = cell(1, nU);
            hasVsm = false;
            for u = 1:nU
                bId = round(TgenAC.bus(rowsU(u)));
                r = VFlexPSession.findBusRecord(recs, bId, {'device'});
                if ~isempty(r) && strcmp(r.model, 'VSM')
                    vsmRecs{u} = r;
                    hasVsm = true;
                end
            end
            if hasVsm
                prm.R_out = zeros(1, nU); prm.L_out = zeros(1, nU);
                prm.Hv    = zeros(1, nU); prm.Dv    = zeros(1, nU);
                prm.Rv    = zeros(1, nU); prm.nq    = zeros(1, nU);
                prm.Tf    = 0.01 * ones(1, nU);
                prm.dynamic_VCC = 'Static';
                for u = 1:nU
                    mb = double(TgenAC.mBase(rowsU(u)));
                    if mb <= 0, mb = Sb; end
                    vr = vsmRecs{u};
                    if ~isempty(vr)
                        p = vr.params;
                        mb = local_pget(p, 'MBASE', mb);
                        sc = Sb / mb;
                        prm.R_out(u) = local_pget(p, 'R_out', 0.01) * sc;
                        prm.L_out(u) = local_pget(p, 'L_out', 0.15) * sc;
                        prm.Hv(u)    = local_pget(p, 'Hv', 5.0) / sc;
                        prm.Dv(u)    = local_pget(p, 'Dv', 0.0) / sc;
                        prm.Rv(u)    = local_pget(p, 'Rv', 0.05) * sc;
                        prm.nq(u)    = local_pget(p, 'nq', 0.05) * sc;
                        prm.Tf(u)    = local_pget(p, 'Tf', 0.01);
                        if strcmpi(local_cget(vr.complexity, 'VCC', 'STATIC'), 'DYNAMIC')
                            prm.dynamic_VCC = 'Dynamic';
                        end
                    else
                        sc = Sb / mb;
                        prm.R_out(u) = 0.01 * sc; prm.L_out(u) = 0.15 * sc;
                        prm.Hv(u) = 5.0 / sc; prm.Dv(u) = 0.0; prm.Rv(u) = 0.05 * sc; prm.nq(u) = 0.05 * sc;
                    end
                end
            end
        end

    case 'GFM'
        if hasDyr
            gfmRecs = cell(1, nU);
            hasGfm = false;
            for u = 1:nU
                bId = round(TgenAC.bus(rowsU(u)));
                r = VFlexPSession.findBusRecord(recs, bId, {'device'});
                if ~isempty(r) && strcmp(r.model, 'GFM')
                    gfmRecs{u} = r;
                    hasGfm = true;
                end
            end
            if hasGfm
                prm.R_out  = 0.01 * ones(1, nU); prm.L_out = 0.1 * ones(1, nU);
                prm.mp     = 0.01 * ones(1, nU); prm.mD    = zeros(1, nU);
                prm.nq     = 0.01 * ones(1, nU); prm.Tf    = 0.01 * ones(1, nU);
                prm.VI_L   = 0.1  * ones(1, nU); prm.VI_R  = 0.01 * ones(1, nU);
                prm.WOF_Tf = zeros(1, nU);
                prm.dynamic_VCC = 'Static';
                for u = 1:nU
                    mb = double(TgenAC.mBase(rowsU(u)));
                    if mb <= 0, mb = Sb; end
                    gr = gfmRecs{u};
                    if ~isempty(gr)
                        p = gr.params;
                        mb = local_pget(p, 'MBASE', mb);
                        sc = Sb / mb;
                        prm.R_out(u)  = local_pget(p, 'R_out', 0.01) * sc;
                        prm.L_out(u)  = local_pget(p, 'L_out', 0.1) * sc;
                        prm.mp(u)     = local_pget(p, 'mp', 0.01) * sc;
                        prm.mD(u)     = local_pget(p, 'mD', 0.0) * sc;
                        prm.nq(u)     = local_pget(p, 'nq', 0.01) * sc;
                        prm.Tf(u)     = local_pget(p, 'Tf', 0.01);
                        prm.VI_L(u)   = local_pget(p, 'VI_L', 0.1) * sc;
                        prm.VI_R(u)   = local_pget(p, 'VI_R', 0.01) * sc;
                        prm.WOF_Tf(u) = local_pget(p, 'WOF_Tf', 0.0);
                        if strcmpi(local_cget(gr.complexity, 'VCC', 'STATIC'), 'DYNAMIC')
                            prm.dynamic_VCC = 'Dynamic';
                        end
                    end
                end
            end
        end

    case 'GFL'
        if hasDyr
            gflRecs = cell(1, nU);
            hasGfl = false;
            for u = 1:nU
                bId = round(TgenAC.bus(rowsU(u)));
                r = VFlexPSession.findBusRecord(recs, bId, {'device'});
                if ~isempty(r) && strcmp(r.model, 'GFL')
                    gflRecs{u} = r;
                    hasGfl = true;
                end
            end
            if hasGfl
                prm.Lf = zeros(1, nU);   prm.Rf = zeros(1, nU);
                prm.Kint = ones(1, nU);  prm.Tint = 0.1 * ones(1, nU); prm.Bint = ones(1, nU);
                prm.kpi = ones(1, nU);   prm.kii = 10.0 * ones(1, nU);
                prm.dynamic_CC = 'Dynamic'; prm.dynamic_PLL = 'Dynamic';
                for u = 1:nU
                    mb = double(TgenAC.mBase(rowsU(u)));
                    if mb <= 0, mb = Sb; end
                    gr = gflRecs{u};
                    if ~isempty(gr)
                        p = gr.params;
                        mb = local_pget(p, 'MBASE', mb);
                        sc = Sb / mb;
                        prm.Lf(u)   = local_pget(p, 'Lf', 0.01) * sc;
                        prm.Rf(u)   = local_pget(p, 'Rf', 0.001) * sc;
                        prm.Kint(u) = local_pget(p, 'Kint', 1.0);
                        prm.Tint(u) = local_pget(p, 'Tint', 0.1);
                        prm.Bint(u) = local_pget(p, 'Bint', 1.0);
                        prm.kpi(u)  = local_pget(p, 'kpi', 1.0);
                        prm.kii(u)  = local_pget(p, 'kii', 10.0);
                        prm.Tf      = local_pget(p, 'Tf', prm.Tf);
                        if strcmpi(local_cget(gr.complexity, 'CC', 'DYNAMIC'), 'STATIC')
                            prm.dynamic_CC = 'Static';
                        end
                        if strcmpi(local_cget(gr.complexity, 'PLL', 'DYNAMIC'), 'STATIC')
                            prm.dynamic_PLL = 'Static';
                        end
                    else
                        sc = Sb / mb;
                        prm.Lf(u) = 0.01 * sc; prm.Rf(u) = 0.001 * sc;
                    end
                end
            end
        end
end

prm = local_apply_overrides(prm, session.unitOverrides, uName);
end

function prm = local_build_dc_unit_prm(uName, rowsD, TgenDC, Sb, session, manifest)
prm = manifest.units.(uName).defaultPrm;
nD  = numel(rowsD);
hasDyr = ~isempty(session.dyrData) && isstruct(session.dyrData) && isfield(session.dyrData, 'records');
if strcmp(uName, 'DC_PES') && hasDyr
    recs = session.dyrData.records;
    pesRecs = cell(1, nD);
    hasPes = false;
    for u = 1:nD
        bId = round(TgenDC.bus(rowsD(u)));
        r = VFlexPSession.findBusRecord(recs, bId, {'device'});
        if ~isempty(r) && strcmp(r.model, 'DC_PES')
            pesRecs{u} = r;
            hasPes = true;
        end
    end
    if hasPes
        prm.Kp_dc = 10.0 * ones(1, nD);
        prm.Ki_dc = 10.0 * ones(1, nD);
        prm.T_pes = 0.25 * ones(1, nD);
        prm.Kd_dc = zeros(1, nD);
        prm.T_wof = 0.01 * ones(1, nD);
        for u = 1:nD
            mb = double(TgenDC.mBase(rowsD(u)));
            if mb <= 0, mb = Sb; end
            pr = pesRecs{u};
            if ~isempty(pr)
                p = pr.params;
                mb = local_pget(p, 'MBASE', mb);
                sc = Sb / mb;
                prm.Kp_dc(u) = local_pget(p, 'Kp_dc', local_pget(p, 'Kp_DC', 10.0)) / sc;
                prm.Ki_dc(u) = local_pget(p, 'Ki_dc', local_pget(p, 'Ki_DC', 10.0)) / sc;
                prm.T_pes(u) = local_pget(p, 'T_pes', local_pget(p, 'T_PES', local_pget(p, 'tau_prim', 0.25)));
                prm.Kd_dc(u) = local_pget(p, 'Kd_dc', 0.0) / sc;
                prm.T_wof(u) = local_pget(p, 'T_wof', 0.01);
            end
        end
    end
end
prm = local_apply_overrides(prm, session.unitOverrides, uName);
end

function prm = local_build_ic_unit_prm(uName, rowsIC, Tintconv, TFb, Sb, session, manifest)
prm = manifest.units.(uName).defaultPrm;
nC  = numel(rowsIC);
hasDyr = ~isempty(session.dyrData) && isstruct(session.dyrData) && isfield(session.dyrData, 'records');
if hasDyr
    recs = session.dyrData.records;
    if strcmp(uName, 'intconv_GFL')
        icRecs = cell(1, nC);
        hasIC = false;
        for u = 1:nC
            r = VFlexPSession.findBranchRecord(recs, Tintconv.fbus(rowsIC(u)), Tintconv.tbus(rowsIC(u)), 'device');
            if ~isempty(r) && strcmp(r.model, 'INTCONV_GFL')
                icRecs{u} = r;
                hasIC = true;
            end
        end
        if hasIC
            prm.Kint = ones(1, nC); prm.Tint = 0.1 * ones(1, nC); prm.Bint = ones(1, nC);
            prm.kpi  = ones(1, nC); prm.kii  = 10.0 * ones(1, nC); prm.Tf = 0.01 * ones(1, nC);
            prm.Kp_dc = 0.5 * ones(1, nC); prm.Ki_dc = 10.0 * ones(1, nC);
            prm.dynamic_CC = 'Dynamic'; prm.dynamic_PLL = 'Dynamic';
            for u = 1:nC
                ir = icRecs{u};
                if ~isempty(ir)
                    p = ir.params;
                    prm.Kint(u)  = local_pget(p, 'Kint', 1.0);
                    prm.Tint(u)  = local_pget(p, 'Tint', 0.1);
                    prm.Bint(u)  = local_pget(p, 'Bint', 1.0);
                    prm.kpi(u)   = local_pget(p, 'kpi', 1.0);
                    prm.kii(u)   = local_pget(p, 'kii', 10.0);
                    prm.Tf(u)    = local_pget(p, 'Tf', 0.01);
                    prm.Kp_dc(u) = local_pget(p, 'Kp_dc', 0.5);
                    prm.Ki_dc(u) = local_pget(p, 'Ki_dc', 10.0);
                    if strcmpi(local_cget(ir.complexity, 'CC', 'DYNAMIC'), 'STATIC')
                        prm.dynamic_CC = 'Static';
                    end
                    if strcmpi(local_cget(ir.complexity, 'PLL', 'DYNAMIC'), 'STATIC')
                        prm.dynamic_PLL = 'Static';
                    end
                end
            end
        end

    elseif strcmp(uName, 'intconv_GFM')
        icRecs = cell(1, nC);
        hasIC = false;
        for u = 1:nC
            r = VFlexPSession.findBranchRecord(recs, Tintconv.fbus(rowsIC(u)), Tintconv.tbus(rowsIC(u)), 'device');
            if ~isempty(r) && strcmp(r.model, 'INTCONV_GFM')
                icRecs{u} = r;
                hasIC = true;
            end
        end
        if hasIC
            prm.R_out  = 0.01 * ones(1, nC); prm.L_out = 0.1 * ones(1, nC);
            prm.mp     = 0.01 * ones(1, nC); prm.nq    = 0.01 * ones(1, nC);
            prm.Tf     = 0.01 * ones(1, nC);
            prm.VI_L   = 0.1  * ones(1, nC); prm.VI_R  = 0.01 * ones(1, nC);
            prm.WOF_Tf = zeros(1, nC);
            prm.dynamic_VCC = 'Static';
            for u = 1:nC
                ir = icRecs{u};
                if ~isempty(ir)
                    p = ir.params;
                    mb = local_pget(p, 'MBASE', Sb);
                    if mb <= 0, mb = Sb; end
                    sc = Sb / mb;
                    prm.R_out(u)  = local_pget(p, 'R_out', 0.01) * sc;
                    prm.L_out(u)  = local_pget(p, 'L_out', 0.1) * sc;
                    prm.mp(u)     = local_pget(p, 'mp', 0.01) * sc;
                    prm.nq(u)     = local_pget(p, 'nq', 0.01) * sc;
                    prm.Tf(u)     = local_pget(p, 'Tf', 0.01);
                    prm.VI_L(u)   = local_pget(p, 'VI_L', 0.1) * sc;
                    prm.VI_R(u)   = local_pget(p, 'VI_R', 0.01) * sc;
                    prm.WOF_Tf(u) = local_pget(p, 'WOF_Tf', 0.0);
                    if strcmpi(local_cget(ir.complexity, 'VCC', 'STATIC'), 'DYNAMIC')
                        prm.dynamic_VCC = 'Dynamic';
                    end
                end
            end
        end

    elseif strcmp(uName, 'intconv_VSM')
        icRecs = cell(1, nC);
        hasIC = false;
        for u = 1:nC
            r = VFlexPSession.findBranchRecord(recs, Tintconv.fbus(rowsIC(u)), Tintconv.tbus(rowsIC(u)), 'device');
            if ~isempty(r) && strcmp(r.model, 'INTCONV_VSM')
                icRecs{u} = r;
                hasIC = true;
            end
        end
        if hasIC
            prm.R_out = reshape(double(Tintconv.r(rowsIC)), 1, []);
            prm.L_out = reshape(double(Tintconv.x(rowsIC)), 1, []);
            prm.Hv = inf(1, nC);         prm.Dv = zeros(1, nC);
            prm.Rv = 0.05 * ones(1, nC); prm.nq = 0.05 * ones(1, nC);
            prm.Tf = 0.01 * ones(1, nC);
            prm.VI_R = zeros(1, nC);     prm.VI_L = zeros(1, nC); prm.WOF_Tf = zeros(1, nC);
            prm.Lf = 0.1 * ones(1, nC);  prm.Rf = 0.01 * ones(1, nC);
            prm.Cfi = 0.1 * ones(1, nC); prm.Rcf = 1e6 * ones(1, nC);
            prm.Kint = ones(1, nC);      prm.Tint = 0.1 * ones(1, nC); prm.Bint = ones(1, nC);
            prm.Kext = ones(1, nC);      prm.Text = 0.1 * ones(1, nC); prm.Bext = ones(1, nC);
            prm.Fi = ones(1, nC);
            prm.Kp_AC = zeros(1, nC);    prm.Kd_AC = zeros(1, nC); prm.N_dc = 100.0 * ones(1, nC);
            prm.Kp_theta = zeros(1, nC); prm.Kp_w  = zeros(1, nC); prm.Ki_w = zeros(1, nC);
            prm.k_p_swing = ones(1, nC);
            prm.dynamic_VCC = 'Static';
            wb = 2 * pi * double(TFb.Fb(1));
            for u = 1:nC
                ir = icRecs{u};
                if ~isempty(ir)
                    p = ir.params;
                    mb = local_pget(p, 'MBASE', Sb);
                    if mb <= 0, mb = Sb; end
                    sc = Sb / mb;
                    prm.R_out(u)  = local_pget(p, 'R_out', prm.R_out(u) / sc) * sc;
                    prm.L_out(u)  = local_pget(p, 'L_out', prm.L_out(u) / sc) * sc;
                    prm.Hv(u)     = local_pget(p, 'Hv', inf) / sc;
                    prm.Dv(u)     = local_pget(p, 'Dv', 0.0) / sc;
                    prm.Rv(u)     = local_pget(p, 'Rv', 0.05) * sc;
                    prm.nq(u)     = local_pget(p, 'nq', 0.05) * sc;
                    prm.Tf(u)     = local_pget(p, 'Tf', 0.01);
                    prm.VI_R(u)   = local_pget(p, 'VI_R', 0.0) * sc;
                    prm.VI_L(u)   = local_pget(p, 'VI_L', 0.0) * sc;
                    prm.WOF_Tf(u) = local_pget(p, 'WOF_Tf', 0.0);
                    prm.Lf(u)     = local_pget(p, 'Lf', 0.1) * sc;
                    prm.Rf(u)     = local_pget(p, 'Rf', 0.01) * sc;
                    prm.Cfi(u)    = local_pget(p, 'Cfi', 0.1) / sc;
                    prm.Rcf(u)    = local_pget(p, 'Rcf', 1e6) * sc;
                    prm.Kint(u)   = local_pget(p, 'Kint', 1.0);
                    prm.Tint(u)   = local_pget(p, 'Tint', 0.1);
                    prm.Bint(u)   = local_pget(p, 'Bint', 1.0);
                    prm.Kext(u)   = local_pget(p, 'Kext', 1.0);
                    prm.Text(u)   = local_pget(p, 'Text', 0.1);
                    prm.Bext(u)   = local_pget(p, 'Bext', 1.0);
                    prm.Fi(u)     = local_pget(p, 'Fi', 1.0);
                    prm.Kp_AC(u) = local_pget(p, 'Kp_AC', local_pget(p, 'Kp_ac', 0.0)) / sc;
                    prm.Kd_AC(u) = local_pget(p, 'Kd_AC', local_pget(p, 'Kd_ac', 0.0)) / sc;
                    prm.N_dc(u)  = local_pget(p, 'N_dc', 100.0);
                    prm.k_p_swing(u) = local_pget(p, 'k_p_swing', 1.0);
                    if isfield(p, 'X_atm')
                        xAtm  = local_pget(p, 'X_atm', 0.0) * sc;
                        zAtm  = local_pget(p, 'zeta_atm', 0.707);
                        wnAtm = local_pget(p, 'wn_atm', 100.0);
                        kpThDef = xAtm;
                        kpWDef  = 4.0 * zAtm * wnAtm * xAtm / wb;
                    else
                        kpThDef = 0.0;
                        kpWDef  = 0.0;
                    end
                    if isfield(p, 'Kp_theta')
                        prm.Kp_theta(u) = local_pget(p, 'Kp_theta', 0.0) * sc;
                    else
                        prm.Kp_theta(u) = kpThDef;
                    end
                    if isfield(p, 'Kp_w')
                        prm.Kp_w(u) = local_pget(p, 'Kp_w', 0.0) * sc;
                    else
                        prm.Kp_w(u) = kpWDef;
                    end
                    prm.Ki_w(u) = local_pget(p, 'Ki_w', 0.0) * sc;
                    if strcmpi(local_cget(ir.complexity, 'VCC', 'STATIC'), 'DYNAMIC')
                        prm.dynamic_VCC = 'Dynamic';
                    end
                end
            end
        end
    end
end
prm = local_apply_overrides(prm, session.unitOverrides, uName);
end

function prm = local_apply_overrides(prm, ovMap, uName)
if isstruct(ovMap) && isfield(ovMap, uName) && isstruct(ovMap.(uName))
    fn = fieldnames(ovMap.(uName));
    for k = 1:numel(fn)
        prm.(fn{k}) = ovMap.(uName).(fn{k});
    end
end
end

function val = local_pget(pstruct, name, defVal)
val = defVal;
if ~isstruct(pstruct), return; end
if isfield(pstruct, name)
    v = pstruct.(name);
    if isnumeric(v) && ~isempty(v), val = double(v(1)); return; end
end
fn = fieldnames(pstruct);
for k = 1:numel(fn)
    if strcmpi(fn{k}, name)
        v = pstruct.(fn{k});
        if isnumeric(v) && ~isempty(v), val = double(v(1)); return; end
    end
end
end

function sval = local_cget(cstruct, name, defStr)
sval = defStr;
if isstruct(cstruct) && isfield(cstruct, name)
    sval = char(string(cstruct.(name)));
end
end

% -------------------------------------------------------------------------
% Analytical & findop Unit Initializers
% -------------------------------------------------------------------------
function op = local_init_op_GFM(op0, prm, vmod, vangle, P, Q)
op = op0;
nU = numel(vmod);
vd_out = real(vmod .* exp(1i * vangle));
vq_out = imag(vmod .* exp(1i * vangle));
id_out = (vd_out .* P + vq_out .* Q) ./ max(vmod.^2, 1e-9);
iq_out = (vq_out .* P - vd_out .* Q) ./ max(vmod.^2, 1e-9);

viR   = local_prm_vec(prm, 'VI_R', 0.01, nU);
viL   = local_prm_vec(prm, 'VI_L', 0.10, nU);
wofTf = local_prm_vec(prm, 'WOF_Tf', 0.0, nU);
rVI   = viR .* double(~(wofTf > 0));
lVI   = viL .* double(~(wofTf > 0));

vc    = (vd_out + 1i * vq_out) + (prm.R_out + 1i * prm.L_out) .* (id_out + 1i * iq_out);
ve    = vc + (rVI + 1i * lVI) .* (id_out + 1i * iq_out);
frame = angle(ve);
e_mag = abs(ve);

vd_c_own   = real(exp(-1i * frame) .* vc);
vq_c_own   = imag(exp(-1i * frame) .* vc);
id_out_own = real(exp(-1i * frame) .* (id_out + 1i * iq_out));
iq_out_own = imag(exp(-1i * frame) .* (id_out + 1i * iq_out));
vd_pcc_own = real(exp(-1i * frame) .* (vd_out + 1i * vq_out));
vq_pcc_own = imag(exp(-1i * frame) .* (vd_out + 1i * vq_out));

id_l = id_out_own + vd_c_own ./ prm.Rcf - prm.Cfi .* vq_c_own;
iq_l = iq_out_own + prm.Cfi .* vd_c_own + vq_c_own ./ prm.Rcf;
PIid = prm.Rf .* id_l ./ prm.Kint - (prm.Bint .* id_l - id_l);
PIiq = prm.Rf .* iq_l ./ prm.Kint - (prm.Bint .* iq_l - iq_l);
PIvd = (id_l - prm.Fi .* id_out_own + prm.Cfi .* vq_c_own) ./ prm.Kext - (prm.Bext .* vd_c_own - vd_c_own);
PIvq = (iq_l - prm.Fi .* iq_out_own - prm.Cfi .* vd_c_own) ./ prm.Kext - (prm.Bext .* vq_c_own - vq_c_own);

op.Pd                         = P;
op.LPF_I                      = vd_pcc_own;
op.LPF1_I                     = vq_pcc_own;
op.LPF_P_I                    = P;
op.LPF_omega_I                = ones(1, nU);
op.LPF_v_I                    = e_mag;
op.PC_LPF_P_Active_Impl_I     = zeros(1, nU);
op.PC_LPF_Q_Active_Impl_I     = zeros(1, nU);
op.PI_synchro_f_I             = zeros(1, nU);
op.PI_synchro_v_I             = zeros(1, nU);
op.RL_out_dynamic_i_d         = id_out;
op.RL_out_dynamic_i_q         = iq_out;
op.VCC_dynamic_RC_dynamic_v_d = vd_c_own;
op.VCC_dynamic_RC_dynamic_v_q = vq_c_own;
op.VCC_dynamic_RL_dynamic_i_d = id_l;
op.VCC_dynamic_RL_dynamic_i_q = iq_l;
op.VCC_dynamic_PIid_I         = PIid;
op.VCC_dynamic_PIiq_I         = PIiq;
op.VCC_dynamic_PIvd_I         = PIvd;
op.VCC_dynamic_PIvq_I         = PIvq;
op.VCC_dynamic_limiting_RC_dynamic_v_d = vd_c_own;
op.VCC_dynamic_limiting_RC_dynamic_v_q = vq_c_own;
op.VCC_dynamic_limiting_RL_dynamic_i_d = id_l;
op.VCC_dynamic_limiting_RL_dynamic_i_q = iq_l;
op.VCC_dynamic_limiting_PIid_I         = PIid;
op.VCC_dynamic_limiting_PIiq_I         = PIiq;
op.VCC_dynamic_limiting_PIvd_I         = PIvd;
op.VCC_dynamic_limiting_PIvq_I         = PIvq;
op.VI_WOF_d_Active_Impl_I              = id_out_own;
op.VI_WOF_q_Active_Impl_I              = iq_out_own;
op.VCC_dynamic_limiting_VI_WOF_d_Active_Impl_I = id_out_own;
op.VCC_dynamic_limiting_VI_WOF_q_Active_Impl_I = iq_out_own;
op.frame                      = frame;
op.iV                         = e_mag - vmod;
op.Integrator                 = 5 * ones(1, nU);
end

function op = local_init_op_GFL(op0, prm, vmod, vangle, P, Q)
op = op0;
nU = numel(vmod);
vd_out = real(vmod .* exp(1i * vangle));
vq_out = imag(vmod .* exp(1i * vangle));
id_g   = (vd_out .* P + vq_out .* Q) ./ max(vmod.^2, 1e-9);
iq_g   = (vq_out .* P - vd_out .* Q) ./ max(vmod.^2, 1e-9);
id_own = real((id_g + 1i * iq_g) .* exp(-1i * vangle));
iq_own = imag((id_g + 1i * iq_g) .* exp(-1i * vangle));

R_out = prm.Rf;
L_out = prm.Lf;
vd_in = vmod + R_out .* id_own - L_out .* iq_own;
vq_in = L_out .* id_own + R_out .* iq_own;
PId   = (vd_in - vmod + L_out .* iq_own) ./ prm.Kint - (prm.Bint .* id_own - id_own);
PIq   = (vq_in - L_out .* id_own) ./ prm.Kint - (prm.Bint .* iq_own - iq_own);

op.CC_dynamic_RL_dynamic_i_d       = id_own;
op.CC_dynamic_RL_dynamic_i_q       = iq_own;
op.CC_dynamic_Subsystem_Integrator  = PId;
op.CC_dynamic_Subsystem1_Integrator = PIq;
op.LPF_I                           = vd_out;
op.LPF1_I                          = vq_out;
op.LPF2_I                          = vmod;
op.LPF3_I                          = zeros(1, nU);
op.LPF_d_Active_Impl_I             = P;
op.LPF_q_Active_Impl_I             = Q;
op.PC_LPF_d_Active_Impl_I          = vmod;
op.PC_LPF_q_Active_Impl_I          = zeros(1, nU);
op.PLL_dynamic_PLL_I               = ones(1, nU);
op.PLL_dynamic_frame               = vangle;
end

function op = local_init_op_intconv_GFM(op0, prm, vmod, vangle, P, Q)
op = op0;
nC = numel(vmod);
vd_out = real(vmod .* exp(1i * vangle));
vq_out = imag(vmod .* exp(1i * vangle));
id_out = (vd_out .* P + vq_out .* Q) ./ max(vmod.^2, 1e-9);
iq_out = (vq_out .* P - vd_out .* Q) ./ max(vmod.^2, 1e-9);

viR   = local_prm_vec(prm, 'VI_R', 0.01, nC);
viL   = local_prm_vec(prm, 'VI_L', 0.10, nC);
wofTf = local_prm_vec(prm, 'WOF_Tf', 0.0, nC);
rVI   = viR .* double(~(wofTf > 0));
lVI   = viL .* double(~(wofTf > 0));

vc    = (vd_out + 1i * vq_out) + (prm.R_out + 1i * prm.L_out) .* (id_out + 1i * iq_out);
ve    = vc + (rVI + 1i * lVI) .* (id_out + 1i * iq_out);
frame = angle(ve);
e_mag = abs(ve);

vd_c_own   = real(exp(-1i * frame) .* vc);
vq_c_own   = imag(exp(-1i * frame) .* vc);
id_out_own = real(exp(-1i * frame) .* (id_out + 1i * iq_out));
iq_out_own = imag(exp(-1i * frame) .* (id_out + 1i * iq_out));
vd_pcc_own = real(exp(-1i * frame) .* (vd_out + 1i * vq_out));
vq_pcc_own = imag(exp(-1i * frame) .* (vd_out + 1i * vq_out));

id_l = id_out_own + vd_c_own ./ prm.Rcf - prm.Cfi .* vq_c_own;
iq_l = iq_out_own + prm.Cfi .* vd_c_own + vq_c_own ./ prm.Rcf;
PIid = prm.Rf .* id_l ./ prm.Kint - (prm.Bint .* id_l - id_l);
PIiq = prm.Rf .* iq_l ./ prm.Kint - (prm.Bint .* iq_l - iq_l);
PIvd = (id_l - prm.Fi .* id_out_own + prm.Cfi .* vq_c_own) ./ prm.Kext - (prm.Bext .* vd_c_own - vd_c_own);
PIvq = (iq_l - prm.Fi .* iq_out_own - prm.Cfi .* vd_c_own) ./ prm.Kext - (prm.Bext .* vq_c_own - vq_c_own);

op.Integrator                 = 5 * ones(1, nC);
op.Integrator1                = e_mag - vmod;
op.LPF_I                      = vd_pcc_own;
op.LPF1_I                     = vq_pcc_own;
op.LPF2_I                     = ones(1, nC);
op.LPF3_I                     = e_mag;
op.PC_LPF_P_Active_Impl_I     = zeros(1, nC);
op.PC_LPF_Q_Active_Impl_I     = zeros(1, nC);
op.PI_I                       = zeros(1, nC);
op.PI1_I                      = zeros(1, nC);
op.RL_out_dynamic_i_d         = id_out;
op.RL_out_dynamic_i_q         = iq_out;
op.VCC_dynamic_RC_dynamic_vd  = vd_c_own;
op.VCC_dynamic_RC_dynamic_vq  = vq_c_own;
op.VCC_dynamic_RL_dynamic_i_d = id_l;
op.VCC_dynamic_RL_dynamic_i_q = iq_l;
op.VCC_dynamic_Subsystem_Integrator  = PIvd;
op.VCC_dynamic_Subsystem1_Integrator = PIvq;
op.VCC_dynamic_Subsystem2_Integrator = PIid;
op.VCC_dynamic_Subsystem3_Integrator = PIiq;
op.VCC_dynamic_limiting_RC_dynamic_vd  = vd_c_own;
op.VCC_dynamic_limiting_RC_dynamic_vq  = vq_c_own;
op.VCC_dynamic_limiting_RL_dynamic_i_d = id_l;
op.VCC_dynamic_limiting_RL_dynamic_i_q = iq_l;
op.VCC_dynamic_limiting_PIv_d_I        = PIvd;
op.VCC_dynamic_limiting_PIv_q_I        = PIvq;
op.VCC_dynamic_limiting_PIi_d_I        = PIid;
op.VCC_dynamic_limiting_PIi_q_I        = PIiq;
op.VI_WOF_d_Active_Impl_I     = id_out_own;
op.VI_WOF_q_Active_Impl_I     = iq_out_own;
op.frame                      = frame;
end

function op = local_init_op_intconv_GFL(op0, prm, dparams, vmod, vangle, P, Q)
op = op0;
nC = numel(vmod);
vd_out = real(vmod .* exp(1i * vangle));
vq_out = imag(vmod .* exp(1i * vangle));
id_g   = (vd_out .* P + vq_out .* Q) ./ max(vmod.^2, 1e-9);
iq_g   = (vq_out .* P - vd_out .* Q) ./ max(vmod.^2, 1e-9);
id_own = real((id_g + 1i * iq_g) .* exp(-1i * vangle));
iq_own = imag((id_g + 1i * iq_g) .* exp(-1i * vangle));

R_out = dparams.Rf;
L_out = dparams.Lf;
vd_in = vmod + R_out .* id_own - L_out .* iq_own;
vq_in = L_out .* id_own + R_out .* iq_own;
PId   = (vd_in - vmod + L_out .* iq_own) ./ prm.Kint - (prm.Bint .* id_own - id_own);
PIq   = (vq_in - L_out .* id_own) ./ prm.Kint - (prm.Bint .* iq_own - iq_own);

op.CC_dynamic_PId_I         = PId;
op.CC_dynamic_PIq_I         = PIq;
op.CC_dynamic_RL_dynamic_i_d = id_own;
op.CC_dynamic_RL_dynamic_i_q = iq_own;
op.PLL_dynamic_frame        = vangle;
op.PLL_dynamic_PLL_I        = ones(1, nC);
op.PC_LPF_d_Active_Impl_I   = vmod;
op.PC_LPF_q_Active_Impl_I   = zeros(1, nC);
op.LPF2_Active_Impl_I       = P;
op.LPF3_Active_Impl_I       = Q;
op.I                        = zeros(1, max(1, nnz(round(dparams.type(:)') == 2)));
end

function op = local_init_op_intconv_VSM(op0, prm, dparams, vmod, vangle, P, Q)
op = op0;
nC = numel(vmod);
vd_out = real(vmod .* exp(1i * vangle));
vq_out = imag(vmod .* exp(1i * vangle));
id_out = (vd_out .* P + vq_out .* Q) ./ max(vmod.^2, 1e-9);
iq_out = (vq_out .* P - vd_out .* Q) ./ max(vmod.^2, 1e-9);

viR   = local_prm_vec(prm, 'VI_R', 0.0, nC);
viL   = local_prm_vec(prm, 'VI_L', 0.0, nC);
wofTf = local_prm_vec(prm, 'WOF_Tf', 0.0, nC);
rVI   = viR .* double(~(wofTf > 0));
lVI   = viL .* double(~(wofTf > 0));

R_out = dparams.Rf;
L_out = dparams.Lf;
vc    = (vd_out + 1i * vq_out) + (R_out + 1i * L_out) .* (id_out + 1i * iq_out);
ve    = vc + (rVI + 1i * lVI) .* (id_out + 1i * iq_out);
frame = angle(ve);
e_mag = abs(ve);
vd_c_own   = real(exp(-1i * frame) .* vc);
vq_c_own   = imag(exp(-1i * frame) .* vc);
id_out_own = real(exp(-1i * frame) .* (id_out + 1i * iq_out));
iq_out_own = imag(exp(-1i * frame) .* (id_out + 1i * iq_out));
vd_pcc_own = real(exp(-1i * frame) .* (vd_out + 1i * vq_out));
vq_pcc_own = imag(exp(-1i * frame) .* (vd_out + 1i * vq_out));

id_l = id_out_own + vd_c_own ./ prm.Rcf - prm.Cfi .* vq_c_own;
iq_l = iq_out_own + prm.Cfi .* vd_c_own + vq_c_own ./ prm.Rcf;
PIid = prm.Rf .* id_l ./ prm.Kint - (prm.Bint .* id_l - id_l);
PIiq = prm.Rf .* iq_l ./ prm.Kint - (prm.Bint .* iq_l - iq_l);
PIvd = (id_l - prm.Fi .* id_out_own + prm.Cfi .* vq_c_own) ./ prm.Kext - (prm.Bext .* vd_c_own - vd_c_own);
PIvq = (iq_l - prm.Fi .* iq_out_own - prm.Cfi .* vd_c_own) ./ prm.Kext - (prm.Bext .* vq_c_own - vq_c_own);

op.RL_out_dynamic_i_d                   = id_out;
op.RL_out_dynamic_i_q                   = iq_out;
op.VCC_dynamic_PIid_I                   = PIid;
op.VCC_dynamic_PIiq_I                   = PIiq;
op.VCC_dynamic_PIvd_Integrator          = PIvd;
op.VCC_dynamic_PIvq_Integrator          = PIvq;
op.VCC_dynamic_RC_dynamic_v_d           = vd_c_own;
op.VCC_dynamic_RC_dynamic_v_q           = vq_c_own;
op.VCC_dynamic_RL_dynamic_i_d           = id_l;
op.VCC_dynamic_RL_dynamic_i_q           = iq_l;
op.VCC_dynamic_limiting_P_PIid_I        = PIid;
op.VCC_dynamic_limiting_P_PIiq_I        = PIiq;
op.VCC_dynamic_limiting_P_PIvd_I        = PIvd;
op.VCC_dynamic_limiting_P_PIvq_I        = PIvq;
op.VCC_dynamic_limiting_P_RC_dynamic_v_d = vd_c_own;
op.VCC_dynamic_limiting_P_RC_dynamic_v_q = vq_c_own;
op.VCC_dynamic_limiting_P_RL_dynamic_i_d = id_l;
op.VCC_dynamic_limiting_P_RL_dynamic_i_q = iq_l;
op.LPF_I                                = vd_pcc_own;
op.LPF1_I                               = vq_pcc_own;
op.VI_WOF_d_Active_Impl_I               = id_out_own;
op.VI_WOF_q_Active_Impl_I               = iq_out_own;
op.frame                                = frame;
op.PC_omega_i                           = ones(1, nC);
op.PC_IP                                = zeros(1, nC);
op.PC_LPF_P_Active_Impl_I               = P;
op.PC_LPF_Q_Active_Impl_I               = Q;
op.PC_LPF_Q1_I                          = Q;
op.LPF3_I                               = P;
op.LPF5_I                               = e_mag;
op.iv                                   = e_mag - vmod;
end

function vec = local_prm_vec(prm, fld, defVal, n)
if isstruct(prm) && isfield(prm, fld) && isnumeric(prm.(fld)) && ~isempty(prm.(fld))
    raw = double(prm.(fld));
    if isscalar(raw)
        vec = raw * ones(1, n);
    else
        vec = reshape(raw, 1, []);
    end
else
    vec = defVal * ones(1, n);
end
end

function [opOut, warnMsg] = local_findop_ac_unit(uName, unitMeta, prmU, op0, v, theta, P, Q, mdlWks, omega_base, Sb)
opOut   = op0;
warnMsg = "";
nU      = numel(v);

if ~bdIsLoaded(unitMeta.libName)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    load_system(fullfile(repoRoot, 'Simulink', 'lib', 'generators', [unitMeta.libName '.slx']));
end

h = new_system();
tmpMdl = get_param(h, 'Name');
cleanTmp = onCleanup(@() bdclose(tmpMdl));
tempWks = get_param(h, 'ModelWorkspace');

% Populate tempWks with clean local routing ([uName '_init_local'] = true)
assignin(tempWks, [uName '_init_local'], true);
assignin(tempWks, [uName 's_buses'], (1:nU)');
assignin(tempWks, [uName '_idx'], mdlWks.getVariable([uName '_idx']));
assignin(tempWks, [uName 's_toM'], mdlWks.getVariable([uName 's_toM']));
assignin(tempWks, ['prm_' uName], prmU);
assignin(tempWks, ['op_' uName], op0);
for varCopy = {'buses_lf','fref','Nbuses','gen_static_params','iPs_lf','iQs_lf','Ps_lf','Qs_lf','From_buses'}
    vn = varCopy{1};
    assignin(tempWks, vn, mdlWks.getVariable(vn));
end
assignin(tempWks, 'Sb', Sb);
assignin(tempWks, 'omega_base', omega_base);

add_block('simulink/Signal Routing/Goto', [tmpMdl '/omega_ref'], 'GotoTag', 'omega_ref', 'TagVisibility', 'global');
add_block('simulink/Signal Routing/Goto', [tmpMdl '/omega_base_goto'], 'GotoTag', 'omega_base', 'TagVisibility', 'global');
add_block('simulink/Signal Routing/Selector', [tmpMdl '/selector1'], ...
    'NumberOfDimensions', '2', 'IndexOptionArray', {'Index vector (dialog)', 'Select all'}, ...
    'IndexParamArray', {'[1:4]', '1'});
add_block('simulink/Signal Routing/Selector', [tmpMdl '/selector2'], ...
    'NumberOfDimensions', '2', 'IndexOptionArray', {'Index vector (dialog)', 'Select all'}, ...
    'IndexParamArray', {'6', '1'});
add_block('simulink/Math Operations/Matrix Concatenate', [tmpMdl '/Matrixconcatenate'], ...
    'NumInputs', '3', 'Mode', 'Multidimensional array', 'ConcatenateDimension', '1');
add_block('simulink/Sources/Constant', [tmpMdl '/one'], ...
    'Value', 'ones(1,length(omega_base))', 'VectorParams1D', 'off');
add_block('simulink/Sources/Constant', [tmpMdl '/omega_base'], ...
    'Value', 'omega_base', 'VectorParams1D', 'off');
add_block('simulink/Sources/Constant', [tmpMdl '/onev'], ...
    'Value', 'ones(1,length(buses_lf))', 'VectorParams1D', 'off');
add_line(tmpMdl, 'one/1', 'omega_ref/1');
add_line(tmpMdl, 'omega_base/1', 'omega_base_goto/1');

uBlk = [tmpMdl '/' uName];
add_block(unitMeta.refBlock, uBlk, 'LinkStatus', 'none');
local_bind_edit_mask_params(uBlk, ['prm_' uName], prmU);

% Set initial guess on active integrators inside temporary unlinked unit
intBlks = find_system(uBlk, 'LookUnderMasks', 'all', ...
    'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
for ii = 1:numel(intBlks)
    if strcmp(uName, 'MGES') && endsWith(intBlks{ii}, '/AVR_LeadLag/x_f') && isfield(prmU, 'T_1')
        nSub = max(1, nnz(prmU.T_1 ~= 0));
        set_param(intBlks{ii}, 'InitialCondition', sprintf('ones(1, %d)', nSub));
    else
        set_param(intBlks{ii}, 'InitialCondition', sprintf('ones(1, %d)', nU));
    end
end

vd = real(v .* exp(1i * theta));
vq = imag(v .* exp(1i * theta));
in_v = [vd; vq; zeros(1, nU)];
assignin(tempWks, 'in_v', in_v);

buses_lf = mdlWks.getVariable('buses_lf');
block_idx = mdlWks.getVariable([uName '_idx']);
elec_in = ones(6, numel(buses_lf));
elec_in(1, block_idx) = 1.0;
elec_in(2, block_idx) = P;
elec_in(3, block_idx) = v;
elec_in(4, block_idx) = Q;
assignin(tempWks, 'elec_in', elec_in);

portsInfo = get_param(uBlk, 'Ports');
nOutOrig  = portsInfo(2);

add_block('simulink/Sources/Constant', [tmpMdl '/v'], 'Value', 'in_v');
add_block('simulink/Sources/Constant', [tmpMdl '/elec'], 'Value', 'elec_in');
add_line(tmpMdl, 'v/1', [uName '/2']);

add_block('simulink/Sinks/Out1', [uBlk '/out_trim1']);
add_line(uBlk, 'trim_output/1', 'out_trim1/1');
add_block('simulink/Sinks/Out1', [uBlk '/out_trim2']);
add_line(uBlk, 'trim_output/2', 'out_trim2/1');
add_block('simulink/Sinks/Out1', [tmpMdl '/out_trim1']);
add_line(tmpMdl, sprintf('%s/%d', uName, nOutOrig + 1), 'out_trim1/1');
add_block('simulink/Sinks/Out1', [tmpMdl '/out_trim2']);
add_line(tmpMdl, sprintf('%s/%d', uName, nOutOrig + 2), 'out_trim2/1');

add_line(tmpMdl, 'elec/1', 'selector1/1');
add_line(tmpMdl, 'elec/1', 'selector2/1');
add_line(tmpMdl, 'selector1/1', 'Matrixconcatenate/1');
add_line(tmpMdl, 'onev/1', 'Matrixconcatenate/2');
add_line(tmpMdl, 'selector2/1', 'Matrixconcatenate/3');
add_line(tmpMdl, 'Matrixconcatenate/1', [uName '/1']);

opspec = operspec(tmpMdl);
opspec.Outputs(1).Known = true(1, nU);
opspec.Outputs(1).y     = (vd .* P + vq .* Q) ./ max(v.^2, 1e-9);
opspec.Outputs(2).Known = true(1, nU);
opspec.Outputs(2).y     = (vq .* P - vd .* Q) ./ max(v.^2, 1e-9);

for s = 1:numel(opspec.States)
    bPath = opspec.States(s).Block;
    szSt  = size(opspec.States(s).x);
    if endsWith(bPath, '/theta') && numel(opspec.States(s).x) == nU
        thGuess = reshape(theta + 0.4, szSt);
        opspec.States(s).x   = thGuess;
        opspec.States(s).Min = reshape(theta - pi/2, szSt);
        opspec.States(s).Max = reshape(theta + 5*pi/6, szSt);
    elseif endsWith(bPath, '/psifd')
        opspec.States(s).x   = ones(szSt);
        opspec.States(s).Min = zeros(szSt);
    elseif endsWith(bPath, '/AVR/Exciter/x_a') || endsWith(bPath, '/AVR_Exciter/x_a')
        opspec.States(s).x   = 1e-3 * ones(szSt);
        opspec.States(s).Min = zeros(szSt);
    elseif endsWith(bPath, '/AVR_LeadLag/x_f')
        opspec.States(s).x   = zeros(szSt);
    elseif (endsWith(bPath, '/psis_q') || endsWith(bPath, '/psikq1') || endsWith(bPath, '/psikq2')) && numel(opspec.States(s).x) == nU
        opspec.States(s).x   = reshape(-v, szSt);
    elseif endsWith(bPath, '/psis_d')
        opspec.States(s).x   = zeros(szSt);
    end
end

optOp = findopOptions( ...
    'OptimizerType', 'graddescent-elim', ...
    'DisplayReport', 'off', ...
    'OptimizationOptions', optimset('MaxFunEvals', 3000 * nU, 'MaxIter', 120, 'TolCon', 1e-10, 'TolX', 1e-11, 'Display', 'off'));
optOp.OptimizationOptions.Algorithm = 'interior-point';

[opFound, opRep] = findop(tmpMdl, opspec, optOp);
termStr = string(opRep.TerminationString);
if ~contains(termStr, "successfully", 'IgnoreCase', true)
    warnMsg = sprintf('%s findop termination: %s', uName, termStr);
end

prefix = [uBlk '/'];
for s = 1:numel(opFound.States)
    bPath = opFound.States(s).Block;
    if startsWith(bPath, prefix)
        rel = extractAfter(bPath, prefix);
        fld = vflexp_state_field(rel);
        opOut.(fld) = reshape(double(opFound.States(s).x), 1, []);
    end
end
end

function [opOut, warnMsg] = local_findop_dc_unit(uName, unitMeta, prmD, op0, DC, mdlWks, Sb)
opOut   = op0;
warnMsg = "";
dcBuses = DC.([uName 's_buses']);
nD      = numel(dcBuses);

if ~bdIsLoaded(unitMeta.libName)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    load_system(fullfile(repoRoot, 'Simulink', 'lib', 'generators', [unitMeta.libName '.slx']));
end

h = new_system();
tmpMdl = get_param(h, 'Name');
cleanTmp = onCleanup(@() bdclose(tmpMdl));
tempWks = get_param(h, 'ModelWorkspace');

DCi = struct();
DCi.v     = DC.v(dcBuses);
DCi.Nbuses = nD;
DCi.([uName 's_Pg'])        = DC.([uName 's_Pg']);
DCi.([uName 's_buses'])     = 1:nD;
DCi.([uName 's_toM'])       = eye(nD);
DCi.([uName 's_agent_idx']) = DC.([uName 's_agent_idx']);

assignin(tempWks, [uName '_init_local'], true);
assignin(tempWks, 'DC', DCi);
assignin(tempWks, 'Sb', Sb);
assignin(tempWks, ['prm_' uName], prmD);
assignin(tempWks, ['op_' uName], op0);

buses_lf = mdlWks.getVariable('buses_lf');
assignin(tempWks, 'buses_lf', buses_lf);
Nag = numel(buses_lf);
ref_in = [ones(1, Nag); ...
          mdlWks.getVariable('Ps_lf'); ...
          mdlWks.getVariable('vs_lf'); ...
          mdlWks.getVariable('Qs_lf'); ...
          ones(1, Nag); ...
          ones(1, Nag)];
assignin(tempWks, 'ref_in', ref_in);

uBlk = [tmpMdl '/' uName];
add_block(unitMeta.refBlock, uBlk, 'LinkStatus', 'none');
local_bind_edit_mask_params(uBlk, ['prm_' uName], prmD);

intBlks = find_system(uBlk, 'LookUnderMasks', 'all', ...
    'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
for ii = 1:numel(intBlks)
    set_param(intBlks{ii}, 'InitialCondition', sprintf('ones(1, %d)', nD));
end

add_block('simulink/Sources/Constant', [tmpMdl '/v'], 'Value', 'DC.v');
add_line(tmpMdl, 'v/1', [uName '/1']);
add_block('simulink/Sources/Constant', [tmpMdl '/ref'], 'Value', 'ref_in');
add_line(tmpMdl, 'ref/1', [uName '/2']);
add_block('simulink/Sinks/Out1', [tmpMdl '/out_trim1']);
add_line(tmpMdl, [uName '/1'], 'out_trim1/1');

opspec = operspec(tmpMdl);
opspec.Outputs(1).Known = true(1, nD);
opspec.Outputs(1).y     = (DC.([uName 's_Pg']) / Sb) ./ max(DCi.v, 1e-6);

optOp = findopOptions( ...
    'OptimizerType', 'graddescent-elim', ...
    'DisplayReport', 'off', ...
    'OptimizationOptions', optimset('MaxFunEvals', 1000 * nD, 'MaxIter', 80, 'TolCon', 1e-7, 'Display', 'off'));

[opFound, opRep] = findop(tmpMdl, opspec, optOp);
termStr = string(opRep.TerminationString);
if ~contains(termStr, "successfully", 'IgnoreCase', true)
    warnMsg = sprintf('%s findop termination: %s', uName, termStr);
end

prefix = [uBlk '/'];
for s = 1:numel(opFound.States)
    bPath = opFound.States(s).Block;
    if startsWith(bPath, prefix)
        rel = extractAfter(bPath, prefix);
        fld = vflexp_state_field(rel);
        opOut.(fld) = reshape(double(opFound.States(s).x), 1, []);
    end
end
end

function local_bind_edit_mask_params(blkPath, prmVarName, prmStruct)
mObj = Simulink.Mask.get(blkPath);
if isempty(mObj), return; end
for k = 1:numel(mObj.Parameters)
    p = mObj.Parameters(k);
    if ismember(p.Type, {'edit', 'promote'}) && isfield(prmStruct, p.Name)
        set_param(blkPath, p.Name, sprintf('%s.%s', prmVarName, p.Name));
    end
end
end
