function [g, dg] = opf_branch_pfsh_fcn_fubm(x, mpc, iPfsh, mpopt)
%OPF_BRANCH_PFSH_FCN_fubm  Evaluates Pf control constraints and Jacobian for elements with active power "from" control using Theta_sh
%   [G, DG] = OPF_BRANCH_FLOW_FCN_AAB(X, OM, IL, MPOPT)
%
%   Element Active Power "from" Control equality constraints for Pf = Pf_set.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X     : optimization vector
%     MPC   : MATPOWER case struct
%     IPFSH : vector of branch indices corresponding to elements with
%             Active power control "Pf = Pf_set".
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (Pf, for all iPfsh) 
%          where the flow will be the active power seen from "from" side
%          The constraint is expressed as (Pf-Pf_set = 0).
%          So far only coded in polar coordinates
%     DG : (optional) equality constraint gradients, column j is
%          gradient of G(j)
%
%   Examples:
%       g = opf_branch_pfsh_fcn_aab(x, mpc, iPfsh, mpopt);
%       [g, dg] = opf_branch_pfsh_fcn_aab(x, mpc, iPfsh, mpopt);
%
%   See also OPF_BRANCH_FLOW_FCN, OPF_BRANCH_FLOW_FCN_AAB, OPF_BRANCH_ZERO_FCN_AAB, .
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
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
[baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);

%% identifier of AC/DC grids
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_branch_pfsh_fcn_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
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
if nPfdp
    branch(iPfdp,SHIFT) = ShAngDp*180/pi;  %AAB- Update the data from Theta_shift Droop to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
end
%% Calculation of admittance matrices
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
Yf=Yf(iPfsh,:); %AAB- Use only the branches within the index iPfsh
Yt=Yt(iPfsh,:); %AAB- Use only the branches within the index iPfsh

%% problem dimensions
nb = length(V);         %% number of buses

%% ----- evaluate constraints -----
if nPfsh > 0 %AAB- If there are Pf control elements
    %% compute branch active power flows of the controlled elements
    Pf = real(V(branch(iPfsh, F_BUS)) .* conj(Yf * V));  %% active power injected at "from" bus (p.u.)
    Pfset = branch(iPfsh, PF)./ mpc.baseMVA; %Setting for active power control in [pu]
    g = [Pf - Pfset];    %% branch active power for Power control (Pf -  Pf_set = 0) (from bus)
else
    g = zeros(0,1);
end

%%----- evaluate partials of constraints -----
if nargout > 1
    if nPfsh > 0
        %% compute partials of Flows w.r.t. V, Beq, Theta Shift and ma
        [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(iPfsh,:), Yf, Yt, V, mpopt.opf.v_cartesian);    %AAB-Obtains the derivatives of the Sf and St w.r.t V   - Yf and Yt are constant here.
        [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1, mpopt.opf.v_cartesian); %% w.r.t. Beqz              %AAB-Obtains the derivatives of the Sf and St w.r.t Beq - V remains constant here because Beq is the only variable
        [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2, mpopt.opf.v_cartesian); %% w.r.t. Beqv              %AAB-Obtains the derivatives of the Sf and St w.r.t Beq - V remains constant here because Beq is the only variable
        [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch(iPfsh,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Theta_sh  %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_sh - V remains constant here because Theta_sh is the only variable 
        [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, mpopt.opf.v_cartesian); %% w.r.t. ma                 %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma/tap   is the only variable
        [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4, mpopt.opf.v_cartesian); %% w.r.t. ma                 %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma/tap   is the only variable
        %[dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3, mpopt.opf.v_cartesian); %% w.r.t. Theta_sh Droop     %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_dp - V remains constant here because Theta_dp is the only variable 

        %% Selecting real part, Pf = real(Sf)
        dPf_dV1  = real(dSf_dV1);
        dPf_dV2  = real(dSf_dV2);
        dPf_dBeqz = real(dSf_dBeqz(iPfsh,:));
        dPf_dBeqv = real(dSf_dBeqv(iPfsh,:));
        dPf_dPfsh = real(dSf_dPfsh);
        dPf_dQtma = real(dSf_dQtma(iPfsh,:));
        dPf_dVtma = real(dSf_dVtma(iPfsh,:)); 
        %dPf_dPfdp = real(dSf_dPfdp(iPfsh,:));
        
        %% construct Jacobian of "from" branch flow inequality constraints
        dg = [ dPf_dV1 dPf_dV2 dPf_dBeqz dPf_dBeqv dPf_dPfsh dPf_dQtma dPf_dVtma];   %% "from" flow limit

    else
        dg = sparse(0, 2*nb+nBeqz+nBeqv+nPfsh+nQtma+nVtma);%<<AAB- No Zero Constrained lines Including FUBM- Original: dh = sparse(0, 2*nb);
    end
end