function d2H = opf_branch_flow_hess_fubm(x, lambda, mpc, il, mpopt)
%OPF_BRANCH_FLOW_HESS_FUBM  Evaluates Hessian of branch flow constraints.
%   D2H = OPF_BRANCH_FLOW_HESS_AAB(X, LAMBDA, OM, IL, MPOPT)
%
%   Hessian evaluation function for AC/DC branch flow constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Kuhn-Tucker multipliers on constrained
%              branch flows
%     MPC : MATPOWER case struct
%     IL : vector of branch indices corresponding to branches with
%          flow limits (all others are assumed to be unconstrained).
%          YF and YT contain only the rows corresponding to IL.
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2H : Hessian of AC branch flow constraints.
%
%   Example:
%       d2H = opf_branch_flow_hess(x, lambda, mpc, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_FCN, OPF_BRANCH_FLOW_HESS.
                                           
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
%% Identifiyng the type of flow limit. "S", "P", "I"
lim_type = upper(mpopt.opf.flow_lim(1));
%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
 
%% identifier of AC/DC grids
%%AAB----------------------------------------------------------------------
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses
%%------------------------------------------------------------------------- 
%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
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
    error('opf_power_balance_hess_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
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
Yf=Yf(il,:);
Yt=Yt(il,:);
%%------------------------------------------------------------------------
%% problem dimensions
nb = length(V);         %% number of buses
nl2 = length(il);       %% number of constrained lines

%%----- evaluate Hessian of flow constraints -----
%% keep dimensions of empty matrices/vectors compatible
%% (required to avoid problems when using Knitro
%%  on cases with all lines unconstrained)
nmu = length(lambda) / 2;
if nmu
    muF = lambda(1:nmu);
    muT = lambda((1:nmu)+nmu);
else    %% keep dimensions of empty matrices/vectors compatible
    muF = zeros(0,1);   %% (required to avoid problems when using Knitro
    muT = zeros(0,1);   %%  on cases with all lines unconstrained)
end
if lim_type == 'I'          %% square of current 
    error('opf_branch_flow_hess_aab: Flow Limit Hessian for "I" type has not been coded yet using FUBM');
    %[dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = dIbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian);
    %d2If_dV2 = @(V, mu)d2Ibr_dV2(Yf, V, mu, mpopt.opf.v_cartesian);
    %d2It_dV2 = @(V, mu)d2Ibr_dV2(Yt, V, mu, mpopt.opf.v_cartesian);
    %[Hf11, Hf12, Hf21, Hf22] = d2Abr_dV2(d2If_dV2, dIf_dV1, dIf_dV2, If, V, muF);
    %[Ht11, Ht12, Ht21, Ht22] = d2Abr_dV2(d2It_dV2, dIt_dV1, dIt_dV2, It, V, muT);
else %<-----AAB-OFSBSR Type == P, 2, S
    f = branch(il, F_BUS);    %% list of "from" buses
    t = branch(il, T_BUS);    %% list of "to" buses
    Cf = sparse(1:nl2, f, ones(nl2, 1), nl2, nb);   %% connection matrix for line & from buses
    Ct = sparse(1:nl2, t, ones(nl2, 1), nl2, nb);   %% connection matrix for line & to buses
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian); %AAB-Obtains the derivatives of the Sf and St w.r.t V   - Yf and Yt are constant here.
    %AAB-------------------------------------------------------------------
    [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1, mpopt.opf.v_cartesian); %% w.r.t. Beq               %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable
    [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2, mpopt.opf.v_cartesian); %% w.r.t. Beq               %AAB-Obtains the derivatives of the Sf and St w.r.t Beq      - V remains constant here because Beq      is the only variable
    [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch(il,:), V, 1, mpopt.opf.v_cartesian);  %% w.r.t. Theta_sh          %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_sh - V remains constant here because Theta_sh is the only variable
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch(il,:), V, 2, mpopt.opf.v_cartesian);  %% w.r.t. ma                %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma       is the only variable
    [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch(il,:), V, 4, mpopt.opf.v_cartesian);  %% w.r.t. ma                %AAB-Obtains the derivatives of the Sf and St w.r.t ma       - V remains constant here because ma       is the only variable
    [dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch(il,:), V, 3, mpopt.opf.v_cartesian);  %% w.r.t. Theta_sh Droop    %AAB-Obtains the derivatives of the Sf and St w.r.t Theta_dp - V remains constant here because Theta_dp is the only variable
    
    %----------------------------------------------------------------------
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. V of Sf 
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, mpopt.opf.v_cartesian); %AAB-Anonymus function for the 2nd derivatives w.r.t. V of St 
    %AAB-------------------------------------------------------------------
    d2Sf_dBeqz2 = @(V, mu)d2Sf_dxBeqz2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqz of Sf 
    d2St_dBeqz2 = @(V, mu)d2St_dxBeqz2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqz of St 
    d2Sf_dBeqv2 = @(V, mu)d2Sf_dxBeqv2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqv of Sf 
    d2St_dBeqv2 = @(V, mu)d2St_dxBeqv2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqv of St 
    d2Sf_dPfsh2 = @(V, mu)d2Sf_dxsh2(branch(il,:), V, mu, mpopt.opf.v_cartesian);     %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of Sf
    d2St_dPfsh2 = @(V, mu)d2St_dxsh2(branch(il,:), V, mu, mpopt.opf.v_cartesian);     %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of St
    d2Sf_dQtma2 = @(V, mu)d2Sf_dxqtma2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. qtma     of Sf
    d2St_dQtma2 = @(V, mu)d2St_dxqtma2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. qtma     of St
    d2Sf_dVtma2 = @(V, mu)d2Sf_dxvtma2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. vtma     of Sf
    d2St_dVtma2 = @(V, mu)d2St_dxvtma2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. vtma     of St
    d2Sf_dPfdp2 = @(V, mu)d2Sf_dxshdp2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of Sf
    d2St_dPfdp2 = @(V, mu)d2St_dxshdp2(branch(il,:), V, mu, mpopt.opf.v_cartesian);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of St
    
    %----------------------------------------------------------------------
    if lim_type == '2'        %% square of real power
         error('opf_branch_flow_hess_aab: Flow Limit Hessian for "P^2" type has not been coded yet using FUBM');
        %[Hf11, Hf12, Hf21, Hf22] = d2Abr_dV2(d2Sf_dV2, real(dSf_dV1), real(dSf_dV2), real(Sf), V, muF);
        %[Ht11, Ht12, Ht21, Ht22] = d2Abr_dV2(d2St_dV2, real(dSt_dV1), real(dSt_dV2), real(St), V, muT);
    elseif lim_type == 'P'    %% real power
        error('opf_branch_flow_hess_aab: Flow Limit Hessian for "P" type has not been coded yet using FUBM');
        %[Hf11, Hf12, Hf21, Hf22] = d2Sf_dV2(V, muF);
        %[Ht11, Ht12, Ht21, Ht22] = d2St_dV2(V, muT);
        %[Hf11, Hf12, Hf21, Hf22] = deal(real(Hf11), real(Hf12), real(Hf21), real(Hf22));
        %[Ht11, Ht12, Ht21, Ht22] = deal(real(Ht11), real(Ht12), real(Ht21), real(Ht22));
    else                      %% square of apparent power
        [Hf11, Hf12, Hf21, Hf22] = d2Abr_dV2(d2Sf_dV2, dSf_dV1, dSf_dV2, Sf, V, muF); 
        [Ht11, Ht12, Ht21, Ht22] = d2Abr_dV2(d2St_dV2, dSt_dV1, dSt_dV2, St, V, muT);
        %AAB---------------------------------------------------------------
        [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Abr_dxBeqz2(d2Sf_dBeqz2, dSf_dV1, dSf_dV2, dSf_dBeqz, Sf, V, muF);
        [Ht13, Ht23, Ht31, Ht32, Ht33] = d2Abr_dxBeqz2(d2St_dBeqz2, dSt_dV1, dSt_dV2, dSt_dBeqz, St, V, muT);
        [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Abr_dxBeqv2(d2Sf_dBeqv2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, Sf, V, muF);
        [Ht14, Ht24, Ht34, Ht41, Ht42, Ht43, Ht44] = d2Abr_dxBeqv2(d2St_dBeqv2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, St, V, muT);
        [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = d2Abr_dxsh2(d2Sf_dPfsh2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, Sf, V, muF); 
        [Ht15, Ht25, Ht35, Ht45, Ht51, Ht52, Ht53, Ht54, Ht55] = d2Abr_dxsh2(d2St_dPfsh2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, St, V, muT); 
        [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = d2Abr_dxqtma2(d2Sf_dQtma2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, Sf, V, muF); 
        [Ht16, Ht26, Ht36, Ht46, Ht56, Ht61, Ht62, Ht63, Ht64, Ht65, Ht66] = d2Abr_dxqtma2(d2St_dQtma2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, St, V, muT); 
        [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = d2Abr_dxvtma2(d2Sf_dVtma2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, dSf_dVtma, Sf, V, muF); 
        [Ht17, Ht27, Ht37, Ht47, Ht57, Ht67, Ht71, Ht72, Ht73, Ht74, Ht75, Ht76, Ht77] = d2Abr_dxvtma2(d2St_dVtma2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, dSt_dVtma, St, V, muT);        
        [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = d2Abr_dxshdp2(d2Sf_dPfdp2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, dSf_dVtma, dSf_dPfdp, Sf, V, muF); 
        [Ht18, Ht28, Ht38, Ht48, Ht58, Ht68, Ht78, Ht81, Ht82, Ht83, Ht84, Ht85, Ht86, Ht87, Ht88] = d2Abr_dxshdp2(d2St_dPfdp2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, dSt_dVtma, dSt_dPfdp, St, V, muT);        
        
        %------------------------------------------------------------------
    end
end
%% construct Hessian
%AAB-----------------------------------------------------------------------
      %Va    Vm    Beqz  Beqv  ShAng Qtma  Vtma  ShAngDp
d2H = [Hf11  Hf12  Hf13  Hf14  Hf15  Hf16  Hf17  Hf18;  
       Hf21  Hf22  Hf23  Hf24  Hf25  Hf26  Hf27  Hf28;
       Hf31  Hf32  Hf33  Hf34  Hf35  Hf36  Hf37  Hf38;
       Hf41  Hf42  Hf43  Hf44  Hf45  Hf46  Hf47  Hf48;
       Hf51  Hf52  Hf53  Hf54  Hf55  Hf56  Hf57  Hf58;
       Hf61  Hf62  Hf63  Hf64  Hf65  Hf66  Hf67  Hf68;
       Hf71  Hf72  Hf73  Hf74  Hf75  Hf76  Hf77  Hf78;
       Hf81  Hf82  Hf83  Hf84  Hf85  Hf86  Hf87  Hf88] + ...
      [Ht11  Ht12  Ht13  Ht14  Ht15  Ht16  Ht17  Ht18;  
       Ht21  Ht22  Ht23  Ht24  Ht25  Ht26  Ht27  Ht28;
       Ht31  Ht32  Ht33  Ht34  Ht35  Ht36  Ht37  Ht38;
       Ht41  Ht42  Ht43  Ht44  Ht45  Ht46  Ht47  Ht48;
       Ht51  Ht52  Ht53  Ht54  Ht55  Ht56  Ht57  Ht58;
       Ht61  Ht62  Ht63  Ht64  Ht65  Ht66  Ht67  Ht68;
       Ht71  Ht72  Ht73  Ht74  Ht75  Ht76  Ht77  Ht78;
       Ht81  Ht82  Ht83  Ht84  Ht85  Ht86  Ht87  Ht88]; %AAB Flow Limit Hessian including FUBM
%--------------------------------------------------------------------------



