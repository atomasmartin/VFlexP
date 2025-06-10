function [g, dg] = opf_branch_pfdp_fcn_fubm(x, mpc, iPfdp, mpopt)
%OPF_BRANCH_PFDP_FCN_fubm  Evaluates Pf control constraints and Jacobian for elements with active power "from" Droop control using Theta_sh
%   [G, DG] = OPF_BRANCH_FLOW_FCN_AAB(X, OM, IL, MPOPT)
%
%   Element Active Power "from" Control equality constraints for Pf = Pf_set.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X     : optimization vector
%     MPC   : MATPOWER case struct
%     IPFDP : vector of branch indices corresponding to elements with
%             Voltage Droop Control "Pf - Pfset = Kdp*(Vmf - Vmfset))".
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
%   Function:
%   br = find(branch(:, BR_STATUS));  %%FUBM- in-service branches
%   brf=branch(:, F_BUS);              %%FUBM- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
%   brt=branch(:, T_BUS);              %%FUBM- to   bus index of all the branches, brt=branch(br, T_BUS); %For in-service branches
%   If= Yf(:, :) * V;                  %%FUBM- complex current injected at "from" bus, Yf(br, :) * V; For in-service branches 
%   It= Yt(:, :) * V;                  %%FUBM- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
%   Sf = V(brf) .* conj(If);           %%FUBM- complex power injected at "from" bus
%   St = V(brt) .* conj(It);           %%FUBM- complex power injected at "to"   bus
%
%   misPfdp = -real(Sf(iPfdp)) + Pfset(iPfdp) + Kdp(iPfdp).* (  Vm(brf(iPfdp)) - Vmfset(iPfdp) )
%
%   Examples:
%       g = opf_branch_pfsh_fcn_aab(x, mpc, iPfdp, mpopt);
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
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_branch_pfsh_fcn_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 6); %AAB- Deals optimisation variables
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
%Yf=Yf(iPfdp,:); %AAB- Use only the branches within the index iPfsh
%Yt=Yt(iPfdp,:); %AAB- Use only the branches within the index iPfsh

%% problem dimensions
nb = length(V);         %% number of buses

%% ----- evaluate constraints -----
if nPfdp > 0 %AAB- If there are Pf control elements
    %% compute branch power flows for the controlled elements
    %br = find(branch(:, BR_STATUS));  %%FUBM- in-service branches
    brf=branch(:, F_BUS);              %%FUBM- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
    brt=branch(:, T_BUS);              %%FUBM- to   bus index of all the branches, brt=branch(br, T_BUS); %For in-service branches
    If= Yf(:, :) * V;                  %%FUBM- complex current injected at "from" bus, Yf(br, :) * V; For in-service branches 
    It= Yt(:, :) * V;                  %%FUBM- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
    Sf = V(brf) .* conj(If);           %%FUBM- complex power injected at "from" bus
    St = V(brt) .* conj(It);           %%FUBM- complex power injected at "to"   bus
    
    %%Save the Voltage-Droop control settings though the branch (   Pf - Pfset = Kdp*(Vmf - Vmfset)  )
    Pfset = branch(:, PF)./ mpc.baseMVA; %Setting for active power control in [pu]
    Kdp   = branch(:,KDP   ); %Voltage Droop Slope   setting for the branch element in p.u.
    Vmfset = branch(:,VF_SET); %Voltage Droop Voltage Setting for the branch element in p.u.

    g = [-real(Sf(iPfdp)) + Pfset(iPfdp) + Kdp(iPfdp).* (  Vm(brf(iPfdp)) - Vmfset(iPfdp) )]; %Droop Pf - Pfset = Kdp*(Vmf - Vmfset)
else
    g = zeros(0,1);
end

%%----- evaluate partials of constraints -----
if nargout > 1
    if nPfdp > 0
        %% compute partials of Flows w.r.t. V, Beq, Theta Shift and ma
        [dPfdp_dV1_all, dPfdp_dV2_all, dPfdp_dPfsh_all, dPfdp_dQfma_all, dPfdp_dBeqz_all,...
            dPfdp_dBeqv_all, dPfdp_dVtma_all, dPfdp_dQtma_all, dPfdp_dPfdp_all] = dPfdp_dx(branch, Yf, Yt, V, 1, mpopt.opf.v_cartesian);
        %% Selecting only the branches that have droop control
        dPfdp_dV1   = dPfdp_dV1_all(iPfdp,:);                %Only Droop control elements iPfdp
        dPfdp_dV2   = dPfdp_dV2_all(iPfdp,:);                %Only Droop control elements iPfdp
        %dPfdp_dPfsh = dPfdp_dPfsh_all(iPfdp,:);              %Only Droop control elements iPfdp
        %dPfdp_dQfma = dPfdp_dQfma_all(iPfdp,:);              %Only Droop control elements iPfdp
        dPfdp_dBeqz = dPfdp_dBeqz_all(iPfdp,:);              %Only Droop control elements iPfdp
        dPfdp_dBeqv = dPfdp_dBeqv_all(iPfdp,:);              %Only Droop control elements iPfdp
        dPfdp_dQtma = dPfdp_dQtma_all(iPfdp,:);              %Only Droop control elements iPfdp
        dPfdp_dVtma = dPfdp_dVtma_all(iPfdp,:);              %Only Droop control elements iPfdp
        dPfdp_dPfdp = dPfdp_dPfdp_all(iPfdp,:);              %Only Droop control elements iPfdp
        
        %% construct Jacobian of "from" branch flow Droop constraints
        dg = [ dPfdp_dV1 dPfdp_dV2 dPfdp_dBeqz dPfdp_dBeqv dPfdp_dQtma dPfdp_dVtma dPfdp_dPfdp];   %% "from" flow limit

    else
        dg = sparse(0, 2*nb+nBeqz+nBeqv+nQtma+nVtma+nPfdp);%<<AAB- No Zero Constrained lines Including FUBM- Original: dh = sparse(0, 2*nb);
    end
end