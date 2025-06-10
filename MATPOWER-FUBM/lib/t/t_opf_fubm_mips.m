function t_opf_fubm_mips(quiet)
%T_OPF_FUBM_MIPS  Tests for MIPS-based ACDC optimal power flow using FUBM formulation.

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 2004-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

%% current mismatch, cartesian V, step-control, linsolver
options = {
    {0, 0, 0, '\'      },
};

num_tests = 56;

t_begin(length(options)*num_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM -extra fields for FUBM
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

casefile = 't_fubm_case_57_14_MTDC_ctrls_opf';
if quiet
    verbose = 0;
else
    verbose = 0;
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

%solver options
mpopt = mpoption('opf.violation', 1e-6);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'mips.gradtol', 1e-8, ...
    'mips.comptol', 1e-8, 'mips.costtol', 1e-9);

for k = 1:length(options)
    if options{k}{1}, bal = 'I';  else, bal = 'S'; end  %% nodal balance
    if options{k}{2}, crd = 'c';  else, crd = 'p'; end  %% V coordinates
    if options{k}{3}, sc = '-sc'; else, sc  = '';  end  %% step control
    t0 = sprintf('MIPS%s (%s,%s,%s) : ', sc, bal, crd, options{k}{4});

    if strcmp(options{k}{4}, 'PARDISO') && ~have_fcn('pardiso')
        t_skip(num_tests, [t0 'PARDISO not available']);
        continue;
    end

    mpopt = mpoption(mpopt, 'opf.current_balance',  options{k}{1}, ...
                            'opf.v_cartesian',      options{k}{2}, ...
                            'mips.step_control',    options{k}{3}, ...
                            'mips.linsolver',       options{k}{4} );

    %% set up indices
    ib_data     = [1:BUS_AREA BASE_KV:VMIN];
    ib_voltage  = [VM VA];
    ib_lam      = [LAM_P LAM_Q];
    ib_mu       = [MU_VMAX MU_VMIN];
    ig_data     = [GEN_BUS QMAX QMIN MBASE:APF];
    ig_disp     = [PG QG VG];
    ig_mu       = (MU_PMAX:MU_QMIN);
    ibr_data    = (1:ANGMAX);
    ibr_flow    = (PF:QT);
    ibr_mu      = [MU_SF MU_ST];
    ibr_angmu   = [MU_ANGMIN MU_ANGMAX];
    ibr_beq     = [BEQ];
    ibr_gsw     = [GSW];
    %% get solved ACDC OPF case from MAT-file
    load soln57-14MTDC_ctrls_fubm_opf;     %% defines bus_soln, gen_soln, branch_soln, f_soln

    %% run OPF
    for s = 0:3
        mpopt = mpoption(mpopt, 'opf.start', s);
        t = sprintf('%s(start=%d): ', t0, s);
        [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt);
        t_ok(success, [t 'success']);
        t_is(f, f_soln, 3, [t 'f']);
        t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ),  3, [t 'bus data']);
        t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
        t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
        t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
        t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ),  3, [t 'gen data']);
        t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
        t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
        t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ),  3, [t 'branch data']); %this includes Theta_sh and ma
        t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
        t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
        t_is(branch(:,ibr_beq   ), branch_soln(:,ibr_beq   ),  3, [t 'branch Beq']);
        t_is(branch(:,ibr_gsw   ), branch_soln(:,ibr_gsw   ),  3, [t 'branch Gsw']);
    end
    mpopt = mpoption(mpopt, 'opf.start', 0);    %% set 'opf.start' back to default
end

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end

t_end;
