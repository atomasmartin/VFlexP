function [g, dg] = opf_power_balance_fcn_fubm(x, mpc, mpopt)
%OPF_POWER_BALANCE_FCN_FUBM  Evaluates AC/DC power balance constraints and their gradients.
%   [G, DG] = OPF_POWER_BALANCE_FCN_AAB(X, OM, MPOPT)
%
%   Computes the Ybus and then computes the active and reactive power balance
%   equality constraints for AC/DC optimal power flow. 
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (active/reactive power balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_power_balance_fcn_aab(x, mpc, mpopt);
%       [g, dg] = opf_power_balance_fcn_aab(x, mpc, mpopt);
%
%   See also OPF_POWER_BALANCE_FCN, OPF_POWER_BALANCE_HESS
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% Identifier of AC/DC grids
%%AAB--------------------------------------------------------------------- 
%%identifier of Zero Constraint VSCs
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1] (Converters and Phase Shifter Transformers, but no VSCIII)
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift (Converters and Phase Shifter Transformers, but no VSCIII)
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_power_balance_fcn_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 1); %AAB- Deals optimisation variables
    V = Vm .* exp(1j * Va);     %% reconstruct V
    %%---------------------------------------------------------------------
end
%% AAB---------------------------------------------------------------------
%%update mpc.branch with FUBM from x
if nBeqz % AC/DC Formulation
    branch(iBeqz,BEQ) = Beqz; %AAB- Update the data from Beqz to the branch matrix
end
if nBeqv
    branch(iBeqv,BEQ) = Beqv; %AAB- Update the data from Beqv to the branch matrix  
end
if nPfsh
    branch(iPfsh,SHIFT) = ShAng*180/pi;  %AAB- Update the data from Theta_shift to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
end
if nQtma
    branch(iQtma,TAP) = maQt;  %AAB- Update the data from ma/tap to the branch matrix.
end
if nVtma
    branch(iVtma,TAP) = maVt;  %AAB- Update the data from ma/tap to the branch matrix.
end
if nPfdp
    branch(iPfdp,SHIFT) = ShAngDp*180/pi;  %AAB- Update the data from Theta_shift Droop to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
end
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %<<AAB-Ybus calculation with updated variables
%%-------------------------------------------------------------------------
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
%% problem dimensions
nb = length(V);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections

%% ----- evaluate constraints -----
%% put Pg, Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% rebuild Sbus
if mpopt.opf.v_cartesian
    Sbus = makeSbus(baseMVA, bus, gen);             %% net injected power in p.u.
else
    Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.
end

%% evaluate complex power balance mismatches
mis = V .* conj(Ybus * V) - Sbus;

%% assemble active and reactive power balance constraints
g = [ real(mis);    %% active power mismatch
      imag(mis) ];  %% reactive power mismatch

%%----- evaluate constraint gradients -----
if nargout > 1
    %% compute partials of injected bus powers
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, mpopt.opf.v_cartesian);   %% w.r.t. V %AAB- Ybus remains constant here because V is the only variable - dSbus_dV1 = dSbus_dVa & dSbus_dV2 = dSbus_dVm for polar
    %%AAB------------------------------------------------------------------
    [dSbus_dBeqz] = dSbus_dBeq(branch, V, 1, mpopt.opf.v_cartesian); %% w.r.t. Beqz      %AAB- V remains constant here because Beq is the only variable
    [dSbus_dBeqv] = dSbus_dBeq(branch, V, 2, mpopt.opf.v_cartesian); %% w.r.t. Beqv      %AAB- V remains constant here because Beq is the only variable
    [dSbus_dPfsh] = dSbus_dsh(branch, V, 1, mpopt.opf.v_cartesian);  %% w.r.t. Theta_sh  %AAB- V remains constant here because Theta_sh is the only variable
    [dSbus_dQtma] = dSbus_dma(branch, V, 2, mpopt.opf.v_cartesian);  %% w.r.t. ma        %AAB- V remains constant here because ma is the only variable    
    [dSbus_dVtma] = dSbus_dma(branch, V, 4, mpopt.opf.v_cartesian);  %% w.r.t. ma        %AAB- V remains constant here because ma is the only variable    
    [dSbus_dPfdp] = dSbus_dsh(branch, V, 3, mpopt.opf.v_cartesian);  %% w.r.t. Theta_sh  %AAB- V remains constant here because Theta_sh is the only variable
    
    %%---------------------------------------------------------------------
    neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);     %% Pbus w.r.t. Pg        %AAB- Also Qbus w.r.t. Qg, It is negative because the generation is negative in the equation: g = mis = V .* conj(Ybus * V) -Sbusg + Sbusd
    if ~mpopt.opf.v_cartesian
        %% adjust for voltage dependent loads
        [dummy, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm);
        dSbus_dV2 = dSbus_dV2 - neg_dSd_dVm;
    end
    dg = [
        real([dSbus_dV1 dSbus_dV2]) neg_Cg sparse(nb, ng) sparse(nb, nBeqz) sparse(nb, nBeqv) real(dSbus_dPfsh) real(dSbus_dQtma) real(dSbus_dVtma) real(dSbus_dPfdp);  %%<<AAB P mismatch w.r.t V1, V2, Pg, Qg, Beqz, Beqv, Theta_sh, maQt, maVt, Theta_shDp  - Original: real([dSbus_dV1 dSbus_dV2]) neg_Cg sparse(nb, ng);
        imag([dSbus_dV1 dSbus_dV2]) sparse(nb, ng) neg_Cg imag(dSbus_dBeqz) imag(dSbus_dBeqv) imag(dSbus_dPfsh) imag(dSbus_dQtma) imag(dSbus_dVtma) imag(dSbus_dPfdp);  %%<<AAB Q mismatch w.r.t V1, V2, Pg, Qg, Beqz, Beqv, Theta_sh, maQt, maVt, Theta_shDp  - Original: imag([dSbus_dV1 dSbus_dV2]) sparse(nb, ng) neg_Cg;
    ];
end