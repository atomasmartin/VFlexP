function d2G = opf_branch_qtma_hess_fubm(x, lambda, mpc, iQtma, mpopt)
%OPF_BRANCH_QTMA_HESS_FUBM  Evaluates Hessian of branch Qt control constraints for elements with reactive power "to" control.
%   D2G = OPF_BRANCH_QTMA_HESS_AAB(X, LAMBDA, MPC, IQTMA, MPOPT)
%
%   Hessian evaluation function for elements with Qt control (Qt - Qtset = 0).
%
%   Inputs:
%     X      : optimization vector
%     LAMBDA : column vector of Kuhn-Tucker multipliers on constrained
%              branch reactive power Qt_set
%     MPC    : MATPOWER case struct
%     IQTMA : vector of branch indices corresponding to elements with
%             Active power control "Qt = Qt_set".
%     MPOPT  : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of Qt control constraints.
%
%   Example:
%       d2G = opf_branch_qtma_hess_aab(x, lambda, mpc, iQtma, mpopt);
%
%   See also OPF_BRANCH_ZERO_HESS, OPF_BRANCH_FLOW_HESS_AAB, OPF_BRANCH_FLOW_HESS, OPF_BRANCH_PFSH_HESS_AAB.
                                           
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
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)& (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
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
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 4); %AAB- Deals optimisation variables
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
Yt=Yt(iQtma,:);

%% problem dimensions
nb = length(V);           %% number of buses
nl = length(branch(:,1)); %% number of lines
%% ----- evaluate Hessian of flow constraints -----
%%keep dimensions of empty matrices/vectors compatible
%%(required to avoid problems when using Knitro
%%on cases with all lines unconstrained)
nmu = length(lambda);
if nmu
    muT = lambda(1:nmu);
else    %% keep dimensions of empty matrices/vectors compatible
    muT = zeros(0,1);   %% required to avoid problems when using Knitro
end
muTaux=sparse(zeros(nl,1)); %% This aux is to select all only the branches that have the Qtma constraint.
muTaux(iQtma)=muT; %% Fill in the location of the constraint the values of mu

    t = branch(iQtma, T_BUS);    %% list of "to" buses
    Ct = sparse(1:nQtma, t, ones(nQtma, 1), nQtma, nb);   %% connection matrix for line & to buses
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. V of St 
    d2St_dBeqz2 = @(V, mu)d2St_dxBeqz2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of St 
    d2St_dBeqv2 = @(V, mu)d2St_dxBeqv2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Beq of St 
    d2St_dsh2 = @(V, mu)d2St_dxsh2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of St 
    d2St_dqtma2 = @(V, mu)d2St_dxqtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. ma       of St 
    %d2St_dvtma2 = @(V, mu)d2St_dxvtma2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. ma       of St 
    d2St_dshdp2 = @(V, mu)d2St_dxshdp2(branch, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of St 

    [Gt11, Gt12, Gt21, Gt22] = d2St_dV2(V, muT);
    [Gt11, Gt12, Gt21, Gt22] = deal(imag(Gt11), imag(Gt12), imag(Gt21), imag(Gt22));
    [Gt13, Gt23, Gt31, Gt32, Gt33] = d2St_dBeqz2(V,muTaux);
    [Gt13, Gt23, Gt31, Gt32, Gt33] = deal(imag(Gt13), imag(Gt23), imag(Gt31), imag(Gt32), imag(Gt33));
    [Gt14, Gt24, Gt34, Gt41, Gt42, Gt43, Gt44] = d2St_dBeqv2(V,muTaux);
    [Gt14, Gt24, Gt34, Gt41, Gt42, Gt43, Gt44] = deal(imag(Gt14), imag(Gt24), imag(Gt34), imag(Gt41), imag(Gt42), imag(Gt43), imag(Gt44));
    [Gt15, Gt25, Gt35, Gt45, Gt51, Gt52, Gt53, Gt54, Gt55] = d2St_dsh2(V,muTaux);
    [Gt15, Gt25, Gt35, Gt45, Gt51, Gt52, Gt53, Gt54, Gt55] = deal(imag(Gt15), imag(Gt25), imag(Gt35), imag(Gt45), imag(Gt51), imag(Gt52), imag(Gt53), imag(Gt54), imag(Gt55));
    [Gt16, Gt26, Gt36, Gt46, Gt56, Gt61, Gt62, Gt63, Gt64, Gt65, Gt66] = d2St_dqtma2(V,muTaux);
    [Gt16, Gt26, Gt36, Gt46, Gt56, Gt61, Gt62, Gt63, Gt64, Gt65, Gt66] = deal(imag(Gt16), imag(Gt26), imag(Gt36), imag(Gt46), imag(Gt56), imag(Gt61), imag(Gt62), imag(Gt63), imag(Gt64), imag(Gt65), imag(Gt66));
    %[Gt17, Gt27, Gt37, Gt47, Gt57, Gt67, Gt71, Gt72, Gt73, Gt74, Gt75, Gt76, Gt77] = d2St_dvtma2(V,muTaux);
    %[Gt17, Gt27, Gt37, Gt47, Gt57, Gt67, Gt71, Gt72, Gt73, Gt74, Gt75, Gt76, Gt77] = deal(imag(Gt17), imag(Gt27), imag(Gt37), imag(Gt47), imag(Gt57), imag(Gt67), imag(Gt71), imag(Gt72), imag(Gt73), imag(Gt74), imag(Gt75), imag(Gt76), imag(Gt77));
    [Gt18, Gt28, Gt38, Gt48, Gt58, Gt68, Gt78, Gt81, Gt82, Gt83, Gt84, Gt85, Gt86, Gt87, Gt88] = d2St_dshdp2(V,muTaux);
    [Gt18, Gt28, Gt38, Gt48, Gt58, Gt68, Gt78, Gt81, Gt82, Gt83, Gt84, Gt85, Gt86, Gt87, Gt88] = deal(imag(Gt18), imag(Gt28), imag(Gt38), imag(Gt48), imag(Gt58), imag(Gt68), imag(Gt78), imag(Gt81), imag(Gt82), imag(Gt83), imag(Gt84), imag(Gt85), imag(Gt86), imag(Gt87), imag(Gt88));
            
    %% construct Hessian
      %Va    Vm    Beqz  Beqv  ShAng Qtma  ShAngDp 
d2G = [Gt11  Gt12  Gt13  Gt14  Gt15  Gt16  Gt18;  
       Gt21  Gt22  Gt23  Gt24  Gt25  Gt26  Gt28;
       Gt31  Gt32  Gt33  Gt34  Gt35  Gt36  Gt38;
       Gt41  Gt42  Gt43  Gt44  Gt45  Gt46  Gt48;
       Gt51  Gt52  Gt53  Gt54  Gt55  Gt56  Gt58;
       Gt61  Gt62  Gt63  Gt64  Gt65  Gt66  Gt68;
       Gt81  Gt82  Gt83  Gt84  Gt85  Gt86  Gt88];%AAB Reactive Flow Qt Constraint Hessian including FUBM




