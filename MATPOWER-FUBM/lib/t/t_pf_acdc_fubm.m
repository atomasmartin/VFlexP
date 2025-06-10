function t_pf_acdc_fubm(quiet)
%T_PF_ACDC_FUBM  Tests for AC/DC power flow solver.

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

ACDC_alg = {'NR'};
ACDC_name = {'Newton (default, power-polar)'};

t_begin(length(ACDC_alg)*12, quiet);

casefile = 't_fubm_case_30_MTDC_ctrls_pf';
if quiet
    verbose = 0;
else
    verbose = 1;
end
if have_fcn('octave')
    if have_fcn('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end
mpopt0 = mpoption('out.all', 0, 'pf.tol', 1e-9, 'verbose', 0);
mpopt0 = mpoption(mpopt0, 'verbose', verbose);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM -extra fields for FUBM
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% -----  AC power flow  -----
%% get solved AC power flow case from MAT-file
load soln30MTDC_ctrls_fubm_pf;      %% defines bus_soln, gen_soln, branch_soln
soln_vg = load('soln30MTDC_ctrls_fubm_pf_vg');

%% run AC PF
for k = 1:length(ACDC_alg)
    t = sprintf('ACDC PF - %s : ', ACDC_name{k});
    mpopt = mpoption(mpopt0, 'pf.alg', ACDC_alg{k});
    [baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
    t_ok(success, [t 'success']);
    t_is(bus, bus_soln, 6, [t 'bus']);
    t_is(gen, gen_soln, 6, [t 'gen']);
    t_is(branch, branch_soln, 6, [t 'branch']);

    r = runpf(casefile, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus, bus_soln, 6, [t 'bus']);
    t_is(r.gen, gen_soln, 6, [t 'gen']);
    t_is(r.branch, branch_soln, 6, [t 'branch']);

    %% check when Vg ~= 1
    t = sprintf('%s - Vg ~= 1 : ', ACDC_alg{k});
    mpopt = mpoption(mpopt, 'verbose', 0);
    mpc = loadcase(casefile);
    mpc.gen(:, VG) = [1.04; 1.025; 1.025; 1.035; 1.03; 1.027];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus, soln_vg.bus_soln, 6, [t 'bus']);
    t_is(r.gen, soln_vg.gen_soln, 6, [t 'gen']);
    t_is(r.branch, soln_vg.branch_soln, 6, [t 'branch']);

end

t_end;

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end
