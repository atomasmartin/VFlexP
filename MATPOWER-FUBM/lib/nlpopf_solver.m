function [results, success, raw] = nlpopf_solver(om, mpopt)
%NLPOPF_SOLVER  Solves AC optimal power flow using MP-Opt-Model.
%
%   [RESULTS, SUCCESS, RAW] = NLPOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln    (deprecated) 2*nb+2*nl - Pmis, Qmis, Sf, St
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .nle    nonlinear equality constraints
%           .nli    nonlinear inequality constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, MIPS.

%   MATPOWER
%   Copyright (c) 2000-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% options
% mips_opt = mpopt.mips;
% mips_opt.verbose = mpopt.verbose;
% if mips_opt.feastol == 0
%     mips_opt.feastol = mpopt.opf.violation;  %% = MPOPT.opf.violation by default
% end
% if ~isfield(mips_opt, 'cost_mult') || isempty(mips_opt.cost_mult)
%     mips_opt.cost_mult = 1e-4;
% end
% opt = struct('mips_opt', mips_opt);

%% unpack data
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nne, nni] = om.get_idx();

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs

%% options
model = om.problem_type();
opt = mpopt2nlpopt(mpopt, model);

%% try to select an interior initial point, unless requested not to
if mpopt.opf.start < 2
    [x0, xmin, xmax] = om.params_var();     %% init var & bounds
    s = 1;                      %% set init point inside bounds by s
    lb = xmin; ub = xmax;
    lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
    ub(xmax ==  Inf) =  1e10;
    x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - s;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + s;                    %% set just above lower bound
    Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);
    Vmax = min(bus(:, VMAX), 1.5);
    Vmin = max(bus(:, VMIN), 0.5);
    Vm = (Vmax + Vmin) / 2;
    if mpopt.opf.v_cartesian
        V = Vm * exp(1j*Varefs(1));
        x0(vv.i1.Vr:vv.iN.Vr) = real(V);
        x0(vv.i1.Vi:vv.iN.Vi) = imag(V);
    else
        x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
        x0(vv.i1.Vm:vv.iN.Vm) = Vm;         %% voltage magnitudes
        if ny > 0
            ipwl = find(gencost(:, MODEL) == PW_LINEAR);
            c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
            x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
        end
    end
    opt.x0 = x0;
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%% Identify FUBM formulation
%the FUBM for DC power flows has not been coded yet
if (size(branch,2) < ALPH3) 
    fubm = 0; %Its not a fubm formulation
else
    fubm = 1; %Its a fubm formulation
end

%% identifier of AC/DC grids
%%AAB--------------------------------------------------------------------- 
if fubm
    iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
    nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
    %%identifier of elements with Vf controlled by Beq
    iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
    nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
    iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
    nVscL = length(iVscL); %AAB- Number of VSC with power losses

    %% Identify if grid has controls
    iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)& (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
    nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
    iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
    nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
    iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
    nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
    iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
    nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
    %%------------------------------------------------------------------------- 
else
    nBeqz = 0;
    nBeqv = 0;
    nPfsh = 0;
    nQtma = 0;
    nVtma = 0;
    nVscL = 0;
    nPfdp = 0;
end
%%-----  run opf  -----
[x, f, eflag, output, lambda] = om.solve(opt);
success = (eflag > 0);

%% update solution data
if mpopt.opf.v_cartesian
    Vi = x(vv.i1.Vi:vv.iN.Vi);
    Vr = x(vv.i1.Vr:vv.iN.Vr);
    V = Vr + 1j*Vi;
    Va = angle(V);
    Vm = abs(V);
else
    Va = x(vv.i1.Va:vv.iN.Va);
    Vm = x(vv.i1.Vm:vv.iN.Vm);
    V = Vm .* exp(1j*Va);
end
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);

%AAB-----------------------------------------------------------------------
if fubm
    if nBeqz %%update mpc.branch with FUBM from x
        Beqz = x(vv.i1.Beqz:vv.iN.Beqz);
        branch(iBeqz,BEQ)=Beqz; %AAB- Update the data from Beqz to the branch matrix
    end
    if nBeqv %%update mpc.branch with FUBM from x
        Beqv = x(vv.i1.Beqv:vv.iN.Beqv);
        branch(iBeqv,BEQ)=Beqv; %AAB- Update the data from Beqz to the branch matrix
    end
    if nPfsh %%update mpc.branch with FUBM from x
        ShAng = x(vv.i1.ShAng:vv.iN.ShAng);
        branch(iPfsh,SHIFT)=ShAng*180/pi; %AAB- Update the data from Theta_sh to the branch matrix (the updated value is changed from radians to degrees)
    end
    if nQtma %%update mpc.branch with FUBM from x
        maQt = x(vv.i1.maQt:vv.iN.maQt);
        branch(iQtma,TAP)=maQt; %AAB- Update the data from ma/tap to the branch matrix 
    end
    if nVtma %%update mpc.branch with FUBM from x
        maVt = x(vv.i1.maVt:vv.iN.maVt);
        branch(iVtma,TAP)=maVt; %AAB- Update the data from ma/tap to the branch matrix 
    end
    if nPfdp %%update mpc.branch with FUBM from x
        ShAngDp = x(vv.i1.ShAngDp:vv.iN.ShAngDp);
        branch(iPfdp,SHIFT)=ShAngDp*180/pi; %AAB- Update the data from Theta_sh to the branch matrix (the updated value is changed from radians to degrees)
    end
end

[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %<<AAB-Ybus calculation with updated variables- Original: makeYbus
%% Standard IEC 62751-2 Ploss Correction for VSC losses
if nVscL
    %%compute branch power flows
    brf=branch(:, F_BUS);              %%AAB- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
    It= Yt(:, :) * V;                  %%AAB- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
    %%compute VSC Power Loss
    PLoss_IEC = branch(iVscL,ALPH3).*((abs(It(iVscL))).^2) + branch(iVscL,ALPH2).*((abs(It(iVscL)))) + branch(iVscL,ALPH1); %%AAB- Standard IEC 62751-2 Ploss Correction for VSC losses 
    branch(iVscL,GSW) = PLoss_IEC./(abs(V(brf(iVscL))).^2);    %%AAB- VSC Gsw Update
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %<<AAB-Ybus calculation with updated variables
end
%--------------------------------------------------------------------------

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = Pg * baseMVA;
gen(:, QG) = Qg * baseMVA;
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch flows
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);    %% build admittance matrices
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);   %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);   %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% line constraint is typically on square of limit
%% so we must fix multipliers
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);
if ~isempty(il)
    if upper(mpopt.opf.flow_lim(1)) == 'P'
        muSf(il) = lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf);
        muSt(il) = lambda.ineqnonlin(nni.i1.St:nni.iN.St);
    else
        muSf(il) = 2 * lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf) .* branch(il, RATE_A) / baseMVA;
        muSt(il) = 2 * lambda.ineqnonlin(nni.i1.St:nni.iN.St) .* branch(il, RATE_A) / baseMVA;
    end
end

%% update Lagrange multipliers
if mpopt.opf.v_cartesian
    if om.userdata.veq
        lam = lambda.eqnonlin(nne.i1.Veq:nne.iN.Veq);
        mu_Vmax = zeros(size(lam));
        mu_Vmin = zeros(size(lam));
        mu_Vmax(lam > 0) =  lam(lam > 0);
        mu_Vmin(lam < 0) = -lam(lam < 0);
        bus(om.userdata.veq, MU_VMAX) = mu_Vmax;
        bus(om.userdata.veq, MU_VMIN) = mu_Vmin;
    end
    bus(om.userdata.viq, MU_VMAX) = lambda.ineqnonlin(nni.i1.Vmax:nni.iN.Vmax);
    bus(om.userdata.viq, MU_VMIN) = lambda.ineqnonlin(nni.i1.Vmin:nni.iN.Vmin);
else
    bus(:, MU_VMAX)  = lambda.upper(vv.i1.Vm:vv.iN.Vm);
    bus(:, MU_VMIN)  = lambda.lower(vv.i1.Vm:vv.iN.Vm);
end
gen(:, MU_PMAX)  = lambda.upper(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = lambda.lower(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = lambda.upper(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = lambda.lower(vv.i1.Qg:vv.iN.Qg) / baseMVA;
if mpopt.opf.current_balance
    %% convert current balance shadow prices to equivalent lamP and lamQ
    %% P + jQ = (Vr + jVi) * (M - jN)
    %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
    %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
    %% lamP = df/dP = df/dM * dM/dP + df/dN + dN/dP
    %% lamQ = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
    VV = V ./ (V .* conj(V));   %% V / Vm^2
    VVr = real(VV);
    VVi = imag(VV);
    lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
    lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
    bus(:, LAM_P) = (VVr.*lamM + VVi.*lamN) / baseMVA;
    bus(:, LAM_Q) = (VVi.*lamM - VVr.*lamN) / baseMVA;
else
    bus(:, LAM_P) = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis) / baseMVA;
    bus(:, LAM_Q) = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis) / baseMVA;
end
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

%% package up results
nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

%% extract multipliers for nonlinear constraints
kl = find(lambda.eqnonlin(1:2*nb) < 0);
ku = find(lambda.eqnonlin(1:2*nb) > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -lambda.eqnonlin(kl);
nl_mu_u(ku) =  lambda.eqnonlin(ku);

if isfield(lambda, 'ineqnonlin')
    lam_nli = lambda.ineqnonlin;
else
    lam_nli = [];
end

mu = struct( ...
  'var', struct('l', lambda.lower, 'u', lambda.upper), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'nle', lambda.eqnonlin, ...
  'nli', lam_nli, ...
  'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', eflag, 'output', output);
