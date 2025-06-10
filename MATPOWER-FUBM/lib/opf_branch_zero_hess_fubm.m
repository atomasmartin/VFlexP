function d2G = opf_branch_zero_hess_fubm(x, lambda, mpc, iBeqz, mpopt)
%OPF_BRANCH_ZERO_HESS_FUBM  Evaluates Hessian of branch Qf zero constraints for VSC.
%   D2G = OPF_BRANCH_ZERO_HESS_AAB(X, LAMBDA, OM, IBEQZ, MPOPT)
%
%   Hessian evaluation function for AC/DC branch Qf zero constraints for VSC.
%
%   Inputs:
%     X      : optimization vector
%     LAMBDA : column vector of Kuhn-Tucker multipliers on constrained
%              branch reactive power "zero constraint"
%     MPC    : MATPOWER case struct
%     IBEQZ  : vector of branch indices corresponding to branches acting as
%              VSC with reactive flow limits (Qf=0).
%     MPOPT  : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of AC/DC branch Qf zero constraints.
%
%   Example:
%       d2G = opf_branch_zero_hess_aab(x, lambda, mpc, iBeq, mpopt);
%
%   See also OPF_BRANCH_ZERO_FCN, OPF_BRANCH_FLOW_HESS_AAB, OPF_BRANCH_FLOW_HESS.
                                           
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
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
nBeqv = 0; %AAB- Vdc control with Beq requires an AC/DC grid.
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

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_branch_zero_hess_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 3); %AAB- Deals optimisation variables
    V = Vm .* exp(1j * Va);     %% reconstruct V
    %%---------------------------------------------------------------------
end
%% AAB---------------------------------------------------------------------
%%update mpc.branch with FUBM from x
if nBeqz % AC/DC Formulation
    branch(iBeqz,BEQ)=Beqz; %AAB- Update the data from Beqz to the branch matrix
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
    branch(iPfdp,SHIFT) = ShAngDp*180/pi;  %AAB- Update the data from Theta_shift to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
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
Yf=Yf(iBeqz,:);

%% problem dimensions
nb = length(V);          %% number of buses
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
muFaux=sparse(zeros(nl,1)); %% This aux is to select all only the branches that have the zero constraint.
muFaux(iBeqz)=muF; %% Fill in the location of the constraint the values of mu

    f = branch(iBeqz, F_BUS);    %% list of "from" buses
    Cf = sparse(1:nBeqz, f, ones(nBeqz, 1), nBeqz, nb);   %% connection matrix for line & from buses
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. V of Sf 
    d2Sf_dBeqz2 = @(V, mu)d2Sf_dxBeqz2(branch(iBeqz,:), V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of Sf 
    %d2Sf_dBeqv2 = @(V, mu)d2Sf_dxBeqv2(branch(iBeqz,:), V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of Sf 
    d2Sf_dsh2 = @(V, mu)d2Sf_dxsh2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of Sf 
    d2Sf_dqtma2 = @(V, mu)d2Sf_dxqtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. ma       of Sf 
    d2Sf_dvtma2 = @(V, mu)d2Sf_dxvtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. ma       of Sf 
    d2Sf_dshdp2 = @(V, mu)d2Sf_dxshdp2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of Sf 

    
    [Gf11, Gf12, Gf21, Gf22] = d2Sf_dV2(V, muF);
    [Gf11, Gf12, Gf21, Gf22] = deal(imag(Gf11), imag(Gf12), imag(Gf21), imag(Gf22));
    [Gf13, Gf23, Gf31, Gf32, Gf33] = d2Sf_dBeqz2(V,muF);
    [Gf13, Gf23, Gf31, Gf32, Gf33] = deal(imag(Gf13), imag(Gf23), imag(Gf31), imag(Gf32), imag(Gf33));
    [Gf15, Gf25, Gf35, Gf45, Gf51, Gf52, Gf53, Gf54, Gf55] = d2Sf_dsh2(V,muFaux);
    [Gf15, Gf25, Gf35, Gf45, Gf51, Gf52, Gf53, Gf54, Gf55] = deal(imag(Gf15), imag(Gf25), imag(Gf35), imag(Gf45), imag(Gf51), imag(Gf52), imag(Gf53), imag(Gf54), imag(Gf55));
    [Gf16, Gf26, Gf36, Gf46, Gf56, Gf61, Gf62, Gf63, Gf64, Gf65, Gf66] = d2Sf_dqtma2(V,muFaux);
    [Gf16, Gf26, Gf36, Gf46, Gf56, Gf61, Gf62, Gf63, Gf64, Gf65, Gf66] = deal(imag(Gf16), imag(Gf26), imag(Gf36), imag(Gf46), imag(Gf56), imag(Gf61), imag(Gf62), imag(Gf63), imag(Gf64), imag(Gf65), imag(Gf66));
    [Gf17, Gf27, Gf37, Gf47, Gf57, Gf67, Gf71, Gf72, Gf73, Gf74, Gf75, Gf76, Gf77] = d2Sf_dvtma2(V,muFaux);
    [Gf17, Gf27, Gf37, Gf47, Gf57, Gf67, Gf71, Gf72, Gf73, Gf74, Gf75, Gf76, Gf77] = deal(imag(Gf17), imag(Gf27), imag(Gf37), imag(Gf47), imag(Gf57), imag(Gf67), imag(Gf71), imag(Gf72), imag(Gf73), imag(Gf74), imag(Gf75), imag(Gf76), imag(Gf77));
    [Gf18, Gf28, Gf38, Gf48, Gf58, Gf68, Gf78, Gf81, Gf82, Gf83, Gf84, Gf85, Gf86, Gf87, Gf88] = d2Sf_dshdp2(V,muFaux);
    [Gf18, Gf28, Gf38, Gf48, Gf58, Gf68, Gf78, Gf81, Gf82, Gf83, Gf84, Gf85, Gf86, Gf87, Gf88] = deal(imag(Gf18), imag(Gf28), imag(Gf38), imag(Gf48), imag(Gf58), imag(Gf68), imag(Gf78), imag(Gf81), imag(Gf82), imag(Gf83), imag(Gf84), imag(Gf85), imag(Gf86), imag(Gf87), imag(Gf88));
    
    %% construct Hessian
      %Va    Vm    Beqz  ShAng Qtma  Vtma  ShAngDp
d2G = [Gf11  Gf12  Gf13  Gf15  Gf16  Gf17  Gf18;  
       Gf21  Gf22  Gf23  Gf25  Gf26  Gf27  Gf28;
       Gf31  Gf32  Gf33  Gf35  Gf36  Gf37  Gf38;
       Gf51  Gf52  Gf53  Gf55  Gf56  Gf57  Gf58;
       Gf61  Gf62  Gf63  Gf65  Gf66  Gf67  Gf68;
       Gf71  Gf72  Gf73  Gf75  Gf76  Gf77  Gf78;
       Gf81  Gf82  Gf83  Gf85  Gf86  Gf87  Gf88];%AAB Reactive Flow zero constraint Hessian including FUBM




