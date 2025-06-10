function d2G = opf_branch_pfsh_hess_fubm(x, lambda, mpc, iPfsh, mpopt)
%OPF_BRANCH_PFSH_HESS_FUBM  Evaluates Hessian of branch Pf control constraints for elements with active power "from" control.
%   D2G = OPF_BRANCH_PFSH_HESS_AAB(X, LAMBDA, OM, IBEQZ, MPOPT)
%
%   Hessian evaluation function for elements with Pf control (Pf - Pfset = 0).
%
%   Inputs:
%     X      : optimization vector
%     LAMBDA : column vector of Kuhn-Tucker multipliers on constrained
%              branch active power Pf_set
%     MPC    : MATPOWER case struct
%     IPFSH : vector of branch indices corresponding to elements with
%             Active power control "Pf = Pf_set".
%     MPOPT  : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of Pf control constraints.
%
%   Example:
%       d2G = opf_branch_pfsh_hess_aab(x, lambda, mpc, iPfsh, mpopt);
%
%   See also OPF_BRANCH_ZERO_HESS, OPF_BRANCH_FLOW_HESS_AAB, OPF_BRANCH_FLOW_HESS.
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% identifier of AC/DC grids
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses

%% Identify if grid has controls
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_branch_pfsh_hess_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 5); %AAB- Deals optimisation variables
    V = Vm .* exp(1j * Va);     %% reconstruct V
    %%---------------------------------------------------------------------
end
%% AAB---------------------------------------------------------------------
%%update mpc.branch with FUBM from x
if nBeqz % AC/DC Formulation
    branch(iBeqz,BEQ)=Beqz; %AAB- Update the data from Beqz to the branch matrix
end
if nBeqv
    branch(iBeqv,BEQ)=Beqv; %AAB- Update the data from Beqv to the branch matrix  
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
%if nPfdp
%    branch(iPfdp,SHIFT) = ShAngDp*180/pi;  %AAB- Update the data from Theta_shift to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
%end
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
Yf=Yf(iPfsh,:);

%% problem dimensions
nb = length(V);           %% number of buses
nl = length(branch(:,1)); %% number of lines
%% ----- evaluate Hessian of flow constraints -----
%%keep dimensions of empty matrices/vectors compatible
%%(required to avoid problems when using Knitro
%%on cases with all lines unconstrained)
nmu = length(lambda);
if nmu
    muF = lambda(1:nmu);
else    %% keep dimensions of empty matrices/vectors compatible
    muF = zeros(0,1);   %% required to avoid problems when using Knitro
end
muFaux=sparse(zeros(nl,1)); %% This aux is to select all only the branches that have the Pfsh constraint.
muFaux(iPfsh)=muF; %% Fill in the location of the constraint the values of mu

    f = branch(iPfsh, F_BUS);    %% list of "from" buses
    Cf = sparse(1:nPfsh, f, ones(nPfsh, 1), nPfsh, nb);   %% connection matrix for line & from buses
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. V of Sf 
    d2Sf_dBeqz2 = @(V, mu)d2Sf_dxBeqz2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of Sf 
    d2Sf_dBeqv2 = @(V, mu)d2Sf_dxBeqv2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of Sf 
    d2Sf_dsh2 = @(V, mu)d2Sf_dxsh2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of Sf 
    d2Sf_dqtma2 = @(V, mu)d2Sf_dxqtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. qtma of Sf 
    d2Sf_dvtma2 = @(V, mu)d2Sf_dxvtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. vtma of Sf 
    %d2Sf_dshdp2 = @(V, mu)d2Sf_dxshdp2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of Sf 

    [Gf11, Gf12, Gf21, Gf22] = d2Sf_dV2(V, muF);
    [Gf11, Gf12, Gf21, Gf22] = deal(real(Gf11), real(Gf12), real(Gf21), real(Gf22));
    [Gf13, Gf23, Gf31, Gf32, Gf33] = d2Sf_dBeqz2(V,muFaux);
    [Gf13, Gf23, Gf31, Gf32, Gf33] = deal(real(Gf13), real(Gf23), real(Gf31), real(Gf32), real(Gf33));
    [Gf14, Gf24, Gf34, Gf41, Gf42, Gf43, Gf44] = d2Sf_dBeqv2(V,muFaux);
    [Gf14, Gf24, Gf34, Gf41, Gf42, Gf43, Gf44] = deal(real(Gf14), real(Gf24), real(Gf34), real(Gf41), real(Gf42), real(Gf43), real(Gf44));
    [Gf15, Gf25, Gf35, Gf45, Gf51, Gf52, Gf53, Gf54, Gf55] = d2Sf_dsh2(V,muFaux);
    [Gf15, Gf25, Gf35, Gf45, Gf51, Gf52, Gf53, Gf54, Gf55] = deal(real(Gf15), real(Gf25), real(Gf35), real(Gf45), real(Gf51), real(Gf52), real(Gf53), real(Gf54), real(Gf55));
    [Gf16, Gf26, Gf36, Gf46, Gf56, Gf61, Gf62, Gf63, Gf64, Gf65, Gf66] = d2Sf_dqtma2(V,muFaux);
    [Gf16, Gf26, Gf36, Gf46, Gf56, Gf61, Gf62, Gf63, Gf64, Gf65, Gf66] = deal(real(Gf16), real(Gf26), real(Gf36), real(Gf46), real(Gf56), real(Gf61), real(Gf62), real(Gf63), real(Gf64), real(Gf65), real(Gf66));
    [Gf17, Gf27, Gf37, Gf47, Gf57, Gf67, Gf71, Gf72, Gf73, Gf74, Gf75, Gf76, Gf77] = d2Sf_dvtma2(V,muFaux);
    [Gf17, Gf27, Gf37, Gf47, Gf57, Gf67, Gf71, Gf72, Gf73, Gf74, Gf75, Gf76, Gf77] = deal(real(Gf17), real(Gf27), real(Gf37), real(Gf47), real(Gf57), real(Gf67), real(Gf71), real(Gf72), real(Gf73), real(Gf74), real(Gf75), real(Gf76), real(Gf77));
    %[Gf18, Gf28, Gf38, Gf48, Gf58, Gf68, Gf78, Gf81, Gf82, Gf83, Gf84, Gf85, Gf86, Gf87, Gf88] = d2Sf_dshdp2(V,muFaux);
    %[Gf18, Gf28, Gf38, Gf48, Gf58, Gf68, Gf78, Gf81, Gf82, Gf83, Gf84, Gf85, Gf86, Gf87, Gf88] = deal(real(Gf18), real(Gf28), real(Gf38), real(Gf48), real(Gf58), real(Gf68), real(Gf78), real(Gf81), real(Gf82), real(Gf83), real(Gf84), real(Gf85), real(Gf86), real(Gf87), real(Gf88));
    
    %% construct Hessian
      %Va    Vm    Beqz  Beqv  ShAng Qtma  Vtma  ShAngDp
d2G = [Gf11  Gf12  Gf13  Gf14  Gf15  Gf16  Gf17;  
       Gf21  Gf22  Gf23  Gf24  Gf25  Gf26  Gf27;
       Gf31  Gf32  Gf33  Gf34  Gf35  Gf36  Gf37;
       Gf41  Gf42  Gf43  Gf44  Gf45  Gf46  Gf47;
       Gf51  Gf52  Gf53  Gf54  Gf55  Gf56  Gf57;
       Gf61  Gf62  Gf63  Gf64  Gf65  Gf66  Gf67;
       Gf71  Gf72  Gf73  Gf74  Gf75  Gf76  Gf77];%AAB Active Flow Pf constraint Hessian including FUBM




