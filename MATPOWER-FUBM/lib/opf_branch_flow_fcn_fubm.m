function [h, dh] = opf_branch_flow_fcn_fubm(x, mpc, il, mpopt)
%OPF_BRANCH_FLOW_FCN_FUBM  Evaluates AC/DC branch flow constraints and Jacobian.
%   [H, DH] = OPF_BRANCH_FLOW_FCN_AAB(X, OM, IL, MPOPT)
%
%   Branch Flow Limit inequality constraints for AC/DC optimal power flow.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     IL : vector of branch indices corresponding to branches with
%          flow limits (all others are assumed to be unconstrained).
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          where the flow will be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (flow^2-limit^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow). So far only apparent
%          power is coded for AC/DC grids
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%
%   Examples:
%       h = opf_branch_flow_fcn_aab(x, mpc, il, mpopt);
%       [h, dh] = opf_branch_flow_fcn_aab(x, mpc, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_HESS, OPF_BRANCH_FLOW_FCN.
                                           
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
lim_type = upper(mpopt.opf.flow_lim(1)); %AAB- Branch flow limit S,P or I
[baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);

%% Identify FUBM formulation
%the FUBM for DC power flows has not been coded yet
if (size(branch,2) < KDP) 
    fubm = 0; %Its not a fubm formulation
else
    fubm = 1; %Its a fubm formulation
end
%% identifier of AC/DC grids
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSCII
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1] (Converters and Phase Shifter Transformers, but no VSCIII)
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift (Converters and Phase Shifter Transformers, but no VSCIII)
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_branch_flow_fcn_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 2); %AAB- Deals optimisation variables
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
Yf=Yf(il,:);
Yt=Yt(il,:);
%%-------------------------------------------------------------------------

%% problem dimensions
nb = length(V);         %% number of buses
nl2 = length(il);       %% number of constrained lines

%% ----- evaluate constraints -----
if nl2 > 0 %AAB- If there are branches with power flow limits
    flow_max = branch(il, RATE_A) / baseMVA; %AAB- Vector with the maximum flow allowed in the branch- size(n12,1)
    if lim_type ~= 'P'      %% typically use square of flow %AAB- For type S and I the limit it's squared
        flow_max = flow_max.^2;
    end
    if lim_type == 'I'      %% current magnitude limit, |I|
        If = Yf * V;
        It = Yt * V;
        h = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
              It .* conj(It) - flow_max ];  %% branch current limits (to bus)
    else % S or P Limit
        %% compute branch power flows of the constrained branches
        Sf = V(branch(il, F_BUS)) .* conj(Yf * V);  %% complex power injected at "from" bus (p.u.)
        St = V(branch(il, T_BUS)) .* conj(Yt * V);  %% complex power injected at "to" bus (p.u.)
        if lim_type == '2'                      %% active power limit, P squared (Pan Wei)
            h = [ real(Sf).^2 - flow_max;       %% branch real power limits (from bus)
                  real(St).^2 - flow_max ];     %% branch real power limits (to bus)
        elseif lim_type == 'P'                  %% active power limit, P
            h = [ real(Sf) - flow_max;          %% branch real power limits (from bus)
                  real(St) - flow_max ];        %% branch real power limits (to bus
        else                                    %% apparent power limit, |S|
            h = [ Sf .* conj(Sf) - flow_max;    %% branch apparent power limits (from bus)
                  St .* conj(St) - flow_max ];  %% branch apparent power limits (to bus)
        end
    end
else
    h = zeros(0,1);
end

%%----- evaluate partials of constraints -----
if nargout > 1
    if nl2 > 0
        %% compute partials of Flows w.r.t. V, Beq, Theta_sh, and ma
        if lim_type == 'I'                      %% current
            %AAB-----------------------------------------------------------
            if (nBeqz || nBeqv || nPfsh || nQtma) %FUBM formulation
                error('opf_branch_flow_fcn_aab: FUBM formulation for limit type "I" has not been coded yet.')
                %[dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian); %AAB-Obtains the derivatives of the Sf and St w.r.t V   - Yf and Yt are constant here.
                %%AAB------------------------------------------------------
                %[dFf_dBeqz, dFt_dBeqz] = dIbr_dBeq(branch(il,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Beqz       %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable
                %[dFf_dBeqv, dFt_dBeqv] = dIbr_dBeq(branch(il,:), V, 2, mpopt.opf.v_cartesian); %% w.r.t. Beqv       %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable                
                %[dFf_dPfsh, dFt_dPfsh] = dIbr_dsh (branch(il,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Theta_sh   %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_sh - V remains constant here because Theta_sh is the only variable                
                %[dFf_dQtma, dFt_dQtma] = dIbr_dma (branch(il,:), V, 2, mpopt.opf.v_cartesian); %% w.r.t. ma         %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma       is the only variable                
                %%---------------------------------------------------------
            else     %AC formulation
                [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian);
            end
            %--------------------------------------------------------------
        else                                    %% power
            %AAB-----------------------------------------------------------
            if fubm %FUBM formulation
                [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian); %AAB-Obtains the derivatives of the Sf and St w.r.t V   - Yf and Yt are constant here.
                %%AAB------------------------------------------------------
                [dFf_dBeqz, dFt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Beqz       %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable
                [dFf_dBeqv, dFt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2, mpopt.opf.v_cartesian); %% w.r.t. Beqv       %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable                
                [dFf_dPfsh, dFt_dPfsh] = dSbr_dsh (branch(il,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Theta_sh   %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_sh - V remains constant here because Theta_sh is the only variable                
                [dFf_dQtma, dFt_dQtma] = dSbr_dma (branch(il,:), V, 2, mpopt.opf.v_cartesian); %% w.r.t. ma         %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma       is the only variable                
                [dFf_dVtma, dFt_dVtma] = dSbr_dma (branch(il,:), V, 4, mpopt.opf.v_cartesian); %% w.r.t. ma         %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma       is the only variable             
                [dFf_dPfdp, dFt_dPfdp] = dSbr_dsh (branch(il,:), V, 3, mpopt.opf.v_cartesian); %% w.r.t. Theta_sh   %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_dp - V remains constant here because Theta_dp is the only variable               
                %%---------------------------------------------------------
            else     %AC formulation
                [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian);
            end
            %--------------------------------------------------------------
        end
        if lim_type == 'P' || lim_type == '2'   %% real part of flow (active power) %AAB- Not coded yet for FUBM
            if (nBeqz || nBeqv || nPfsh || nQtma|| nVtma) %FUBM formulation
                dFf_dV1 = real(dFf_dV1);
                dFf_dV2 = real(dFf_dV2);
                dFt_dV1 = real(dFt_dV1);
                dFt_dV2 = real(dFt_dV2);
                %AAB-------------------------------------------------------
                dFf_dBeqz = real(dFf_dBeqz);
                dFf_dBeqv = real(dFf_dBeqv);
                dFf_dPfsh = real(dFf_dPfsh);
                dFf_dQtma = real(dFf_dQtma);
                dFf_dVtma = real(dFf_dVtma);
                dFf_dPfdp = real(dFf_dPfdp);
                
                dFt_dBeqz = real(dFt_dBeqz);
                dFt_dBeqv = real(dFt_dBeqv);
                dFt_dPfsh = real(dFt_dPfsh);
                dFt_dQtma = real(dFt_dQtma);
                dFt_dVtma = real(dFt_dVtma);
                dFt_dPfdp = real(dFt_dPfdp);
                %----------------------------------------------------------
                Ff = real(Ff);
                Ft = real(Ft);
            else %AC formulation
                dFf_dV1 = real(dFf_dV1);
                dFf_dV2 = real(dFf_dV2);
                dFt_dV1 = real(dFt_dV1);
                dFt_dV2 = real(dFt_dV2);
                Ff = real(Ff);
                Ft = real(Ft);
            end
        end
        
        %% Derivatives for "Limit^2 - Flow^2" (For "I" or "S") or "Limit - Flow" (For "P")
        if lim_type == 'P'  %AAB- "Limit - Flow"
            %% active power
            if fubm %FUBM formulation
                [df_dV1, df_dV2, dt_dV1, dt_dV2] = deal(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2);
                [df_dBeqz, df_dBeqv, df_dPfsh, df_dQtma, df_dVtma, df_dPfdp] = deal(dFf_dBeqz, dFf_dBeqv, dFf_dPfsh, dFf_dQtma, dFf_dVtma, dFf_dPfdp);
                [dt_dBeqz, dt_dBeqv, dt_dPfsh, dt_dQtma, dt_dVtma, dt_dPfdp] = deal(dFt_dBeqz, dFt_dBeqv, dFt_dPfsh, dFt_dQtma, dFt_dVtma, dFt_dPfdp);
            else %AC Formulation 
                [df_dV1, df_dV2, dt_dV1, dt_dV2] = deal(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2); 
            end
        else %<-----AAB- "Limit^2 - Flow^2"- OFSBSR
            %AAB-----------------------------------------------------------
            if fubm  %%<-----AAB- FUBM formulation -OFSBSR
                %%squared magnitude of flow (of complex power or current, or real power)
                [df_dV1, df_dV2, dt_dV1, dt_dV2] = ...
                    dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
                [df_dBeqz, dt_dBeqz] = ...
                    dAbr_dBeq(dFf_dBeqz, dFt_dBeqz, Ff, Ft);
                [df_dBeqv, dt_dBeqv] = ...
                        dAbr_dBeq(dFf_dBeqv, dFt_dBeqv, Ff, Ft);   
                [df_dPfsh, dt_dPfsh] = ...
                        dAbr_dsh(dFf_dPfsh, dFt_dPfsh, Ff, Ft); 
                [df_dQtma, dt_dQtma] = ...
                        dAbr_dma(dFf_dQtma, dFt_dQtma, Ff, Ft);  
                [df_dVtma, dt_dVtma] = ...
                        dAbr_dma(dFf_dVtma, dFt_dVtma, Ff, Ft);
                [df_dPfdp, dt_dPfdp] = ...
                        dAbr_dsh(dFf_dPfdp, dFt_dPfdp, Ff, Ft); 
            else     %AC formulation
                %%squared magnitude of flow (of complex power or current, or real power)
                [df_dV1, df_dV2, dt_dV1, dt_dV2] = ...
                    dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
            end
            %--------------------------------------------------------------
        end

        %% construct Jacobian of "from" and "to" branch flow inequality constraints
        %AAB---------------------------------------------------------------
        if fubm %FUBM formulation
            dh = [ df_dV1 df_dV2 df_dBeqz df_dBeqv df_dPfsh df_dQtma df_dVtma df_dPfdp;      %% "from" flow limit
                   dt_dV1 dt_dV2 dt_dBeqz dt_dBeqv dt_dPfsh dt_dQtma dt_dVtma dt_dPfdp];     %%  "to"  flow limit
        else     %AC formulation
            dh = [ df_dV1 df_dV2;                   %% "from" flow limit
                   dt_dV1 dt_dV2 ];                 %%  "to"  flow limit
        end
        %------------------------------------------------------------------
    else
        dh = sparse(0, 2*nb+nBeqz+nBeqv+nPfsh+nQtma+nVtma+nPfdp);%<<AAB- No Constrained lines Including FUBM- Original: dh = sparse(0, 2*nb);
    end
end
