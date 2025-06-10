function [num_Hf18, num_Hf28, num_Hf38, num_Hf48, num_Hf58, num_Hf68, num_Hf78, num_Hf81, num_Hf82, num_Hf83, num_Hf84, num_Hf85, num_Hf86, num_Hf87, num_Hf88,...
          num_Ht18, num_Ht28, num_Ht38, num_Ht48, num_Ht58, num_Ht68, num_Ht78, num_Ht81, num_Ht82, num_Ht83, num_Ht84, num_Ht85, num_Ht86, num_Ht87, num_Ht88] = d2Abr_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXSHDP2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. ShdpVa, ShdpVm, ShdpBeqz, ShdpBeqv, ShdpSh, Shdpqtma, VaShdp, VmShdp, BeqzShdp, BeqvShdp, ShShdp, qtmaShdp, vtmaShdp, ShdpShdp(Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%  [NUM_HF18, NUM_HF28, NUM_HF38, NUM_HF48, NUM_HF58, NUM_HF68, NUM_HF78, NUM_HF81, NUM_HF82, NUM_HF83, NUM_HF84, NUM_HF85, NUM_HF86, NUM_HF87, NUM_HF88,...
%   NUM_HT18, NUM_HT28, NUM_HT38, NUM_HT48, NUM_HT58, NUM_HT68, NUM_HT78, NUM_HT81, NUM_HT82, NUM_HT83, NUM_HT84, NUM_HT85, NUM_HT86, NUM_HT87, NUM_HT88] = D2ABR_DXSHDP2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 30 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H18 = HShdpVa   = d2Abr_dPfdpVa     = d/dPfdp ( (dAbr/dVa    )'*mu )
%   H28 = HShdpVm   = d2Abr_dPfdpVm     = d/dPfdp ( (dAbr/dVm    )'*mu )
%   H38 = HShdpBz   = d2Abr_dPfdpBeqz   = d/dPfdp ( (dAbr/dBeqz  )'*mu )
%   H48 = HShdpBv   = d2Abr_dPfdpBeqv   = d/dPfdp ( (dAbr/dBeqv  )'*mu )
%   H58 = HShdpSh   = d2Abr_dPfdpSh     = d/dPfdp ( (dAbr/dSh    )'*mu )
%   H68 = HShdpSh   = d2Abr_dPfdpqtma   = d/dPfdp ( (dAbr/dqtma  )'*mu )
%   H78 = HShdpSh   = d2Abr_dPfdpvtma   = d/dPfdp ( (dAbr/dvtma  )'*mu )
%   H81 = HVaShdp   = d2Abr_dVaPfdp     = d/dVa   ( (dAbr/dPfdp  )'*mu )
%   H82 = HVmShdp   = d2Abr_dVmPfdp     = d/dVm   ( (dAbr/dPfdp  )'*mu )
%   H83 = HBzShdp   = d2Abr_dBeqzPfdp   = d/dBeqz ( (dAbr/dPfdp  )'*mu )
%   H84 = HBvShdp   = d2Abr_dBeqvPfdp   = d/dBeqv ( (dAbr/dPfdp  )'*mu )
%   H85 = HShShdp   = d2Abr_dShPfdp     = d/dsh   ( (dAbr/dPfdp  )'*mu )
%   H86 = HqtmaShdp = d2Abr_dqtmaPfdp   = d/dqtma ( (dAbr/dPfdp  )'*mu )
%   H87 = HqtmaShdp = d2Abr_dqtmaPfdp   = d/dvtma ( (dAbr/dPfdp  )'*mu )
%   H88 = HvtmaShdp = d2Abr_dvtmaPfdp   = d/dPfdp ( (dAbr/dPfdp  )'*mu )
%   
%   H18 = HShdpVa   = d2Abr_dPfdpVa     =  (  2*real(d2Sbr_dPfdpVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dPfdp   ))'*mu 
%   H28 = HShdpVm   = d2Abr_dPfdpVm     =  (  2*real(d2Sbr_dPfdpVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dPfdp   ))'*mu
%   H38 = HShdpBz   = d2Abr_dPfdpBeqz   =  (  2*real(d2Sbr_dPfdpBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dPfdp   ))'*mu
%   H48 = HShdpBv   = d2Abr_dPfdpBeqv   =  (  2*real(d2Sbr_dPfdpBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dPfdp   ))'*mu
%   H58 = HShdpSh   = d2Abr_dPfdpSh     =  (  2*real(d2Sbr_dPfdpSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dPfdp   ))'*mu
%   H68 = HShdpqtma = d2Abr_dPfdpqtma   =  (  2*real(d2Sbr_dPfdpqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dPfdp   ))'*mu
%   H78 = HShdpqtma = d2Abr_dPfdpqtma   =  (  2*real(d2Sbr_dPfdpvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dPfdp   ))'*mu

%   H81 = HVaShdp   = d2Abr_dVaPfdp     =  (  2*real(d2Sbr_dVaPfdp  )*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dVa     ))'*mu
%   H82 = HVmShdp   = d2Abr_dVmPfdp     =  (  2*real(d2Sbr_dVmPfdp  )*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dVm     ))'*mu
%   H83 = HBzShdp   = d2Abr_dBeqzPfdp   =  (  2*real(d2Sbr_dBeqzPfdp)*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dBeqz   ))'*mu
%   H84 = HBvShdp   = d2Abr_dBeqvPfdp   =  (  2*real(d2Sbr_dBeqvPfdp)*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dBeqv   ))'*mu
%   H85 = HShShdp   = d2Abr_dShPfdp     =  (  2*real(d2Sbr_dShPfdp  )*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dSh     ))'*mu 
%   H86 = HqtmaShdp = d2Abr_dqtmaPfdp   =  (  2*real(d2Sbr_dqtmaPfdp)*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dqtma   ))'*mu 
%   H87 = HvtmaShdp = d2Abr_dvtmaPfdp   =  (  2*real(d2Sbr_dvtmaPfdp)*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dvtma   ))'*mu 
%   H88 = HShdpShdp = d2Abr_dPfdpPfdp   =  (  2*real(d2Sbr_dPfdpPfdp)*conj(Sbr)  +  dSbr_dPfdp *conj(dSbr_dPfdp   ))'*mu 
%
%   Example:
%  [NUM_HF18, NUM_HF28, NUM_HF38, NUM_HF48, NUM_HF58, NUM_HF68, NUM_HF78, NUM_HF81, NUM_HF82, NUM_HF83, NUM_HF84, NUM_HF85, NUM_HF86, NUM_HF87, NUM_HF88,...
%   NUM_HT18, NUM_HT28, NUM_HT38, NUM_HT48, NUM_HT58, NUM_HT68, NUM_HT78, NUM_HT81, NUM_HT82, NUM_HT83, NUM_HT84, NUM_HT85, NUM_HT86, NUM_HT87, NUM_HT88] = D2ABR_DXSHDP2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXSH2.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TNX]  A. Alvarez-Bustos, "AC/DC FUBM and their Derivatives using FUBM
%          Complex Matrix Notation" MATPOWER Technical Note x, Month 20XX.
%             http://www.pserc.cornell.edu/matpower/
%                                           TNX-OPF-Derivatives-FUBM.pdf   %%AAB- Technical note to be written
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM

%% default input args
if nargin < 7
    vcart = 0;      %% default to polar coordinates
end

%% constants
nl = length(branch(:,1)); %% number of lines
nb = length(V);  %% number of buses
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
%Changing Perturbation to Degrees for Theta_sh since inside makeYbus it gets changed to radians.
pertDeg = (pert*180)/pi;
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift

%% Calculation of derivatives
if vcart
    error('d2Abr_dxShdp2Pert: Derivatives of Flow Limit equations w.r.t Theta_dp using Finite Differences Method in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
       
    %Make Ybus, Yf, Yt
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
    %Sbr 1st Derivatives 
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1, vcart);
    [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2, vcart);
    [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1, vcart);
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, vcart);
    [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4, vcart);    
    [dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3, vcart);
    
    %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St); 
    [dAf_dPfsh, dAt_dPfsh] = dAbr_dsh(dSf_dPfsh, dSt_dPfsh, Sf, St); 
    [dAf_dQtma, dAt_dQtma] = dAbr_dma(dSf_dQtma, dSt_dQtma, Sf, St); 
    [dAf_dVtma, dAt_dVtma] = dAbr_dma(dSf_dVtma, dSt_dVtma, Sf, St); 
    [dAf_dPfdp, dAt_dPfdp] = dAbr_dsh(dSf_dPfdp, dSt_dPfdp, Sf, St);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
      
    %Selector of active Beq 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
        
    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
 
    %Selector of active Qt ma/tap 
    QtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    QtmaSel(iQtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Qt
    diagQtmaSel = sparse( diag(QtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Qtma, size [nl,nl]
    
    %Selector of active Vt ma/tap 
    VtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    VtmaSel(iVtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Vt
    diagVtmaSel = sparse( diag(VtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Vtma, size [nl,nl]
    
    %Selector of active Theta_dp
    ShdpSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    ShdpSel(iPfdp) = 1;                       %AAB- Fill the selector with 1 where Theta_dp is active controlling Droop
    diagShdpSel = sparse( diag(ShdpSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShdpShdp, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)  
    d2Af_dPfdpVa   = sparse( zeros(nb,   nPfdp) ); 
    d2Af_dPfdpVm   = sparse( zeros(nb,   nPfdp) );
    d2Af_dPfdpBeqz = sparse( zeros(nBeqz,nPfdp) );
    d2Af_dPfdpBeqv = sparse( zeros(nBeqv,nPfdp) );
    d2Af_dPfdpPfsh = sparse( zeros(nPfsh,nPfdp) );   
    d2Af_dPfdpQtma = sparse( zeros(nQtma,nPfdp) );
    d2Af_dPfdpVtma = sparse( zeros(nVtma,nPfdp) );    
    
    d2Af_dVaPfdp   = sparse( zeros(nPfdp,nb   ) ); 
    d2Af_dVmPfdp   = sparse( zeros(nPfdp,nb   ) );
    d2Af_dBeqzPfdp = sparse( zeros(nPfdp,nBeqz) );
    d2Af_dBeqvPfdp = sparse( zeros(nPfdp,nBeqv) );
    d2Af_dPfshPfdp = sparse( zeros(nPfdp,nPfsh) );   
    d2Af_dQtmaPfdp = sparse( zeros(nPfdp,nQtma) );
    d2Af_dVtmaPfdp = sparse( zeros(nPfdp,nVtma) );    
    
    d2Af_dPfdp2    = sparse( zeros(nPfdp,nPfdp) );

    d2At_dPfdpVa   = sparse( zeros(nb,   nPfdp) ); 
    d2At_dPfdpVm   = sparse( zeros(nb,   nPfdp) );
    d2At_dPfdpBeqz = sparse( zeros(nBeqz,nPfdp) );
    d2At_dPfdpBeqv = sparse( zeros(nBeqv,nPfdp) );
    d2At_dPfdpPfsh = sparse( zeros(nPfsh,nPfdp) );   
    d2At_dPfdpQtma = sparse( zeros(nQtma,nPfdp) );  
    d2At_dPfdpVtma = sparse( zeros(nVtma,nPfdp) );      
    
    d2At_dVaPfdp   = sparse( zeros(nPfdp,nb   ) ); 
    d2At_dVmPfdp   = sparse( zeros(nPfdp,nb   ) );
    d2At_dBeqzPfdp = sparse( zeros(nPfdp,nBeqz) );
    d2At_dBeqvPfdp = sparse( zeros(nPfdp,nBeqv) );
    d2At_dPfshPfdp = sparse( zeros(nPfdp,nPfsh) );   
    d2At_dQtmaPfdp = sparse( zeros(nPfdp,nQtma) );
    d2At_dVtmaPfdp = sparse( zeros(nPfdp,nVtma) );    
   
    d2At_dPfdp2    = sparse( zeros(nPfdp,nPfdp) );
    
    %PfdpVa num_G81
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dPfdp_PertVa, dSt_dPfdp_PertVa] = dSbr_dsh(branch, V1p, 3, vcart); %dSbr_dPfdpPertVa
        [dAf_dPfdp_PertVa, dAt_dPfdp_PertVa] = dAbr_dsh(dSf_dPfdp_PertVa, dSt_dPfdp_PertVa, Sf_PertVa, St_PertVa); %dAbr_dPfdpPertVa       
        d2Af_dVaPfdp(:, k) = (dAf_dPfdp_PertVa - dAf_dPfdp).' * lam / pert; %PfdpVa from, size of [nPfdp, nb]
        d2At_dVaPfdp(:, k) = (dAt_dPfdp_PertVa - dAt_dPfdp).' * lam / pert; %PfdpVa  to , size of [nPfdp, nb]
    end
    %PfdpVm num_G82
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dPfdp_PertVm, dSt_dPfdp_PertVm] = dSbr_dsh(branch, V2p, 3, vcart); %dSbr_dPfdpPertVm
        [dAf_dPfdp_PertVm, dAt_dPfdp_PertVm] = dAbr_dsh(dSf_dPfdp_PertVm, dSt_dPfdp_PertVm, Sf_PertVm, St_PertVm); %dAbr_dPfdpPertVm       
        d2Af_dVmPfdp(:, k) = (dAf_dPfdp_PertVm - dAf_dPfdp).' * lam / pert; %PfdpVm from, size of [nPfdp, nb]
        d2At_dVmPfdp(:, k) = (dAt_dPfdp_PertVm - dAt_dPfdp).' * lam / pert; %PfdpVm  to , size of [nPfdp, nb]
    end
    %PfdpBeqz num_G83
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertBeqz, St_PertBeqz] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertBeqz evaluated in x+pert
        [dSf_dPfdp_PertBeqz, dSt_dPfdp_PertBeqz] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertBeqz
        %dAbr_dPfdpPertBeqz evaluated in x+pert        
        [dAf_dPfdp_PertBeqz, dAt_dPfdp_PertBeqz] = dAbr_dsh(dSf_dPfdp_PertBeqz, dSt_dPfdp_PertBeqz, Sf_PertBeqz, St_PertBeqz); %dAbr_dPfdpPertBeqz              
        %2nd Derivatives of Abr w.r.t. PfdpBeqz
        d2Af_dBeqzPfdp(:, k) = (dAf_dPfdp_PertBeqz - dAf_dPfdp).' * lam / pert; %PfdpBeqz From
        d2At_dBeqzPfdp(:, k) = (dAt_dPfdp_PertBeqz - dAt_dPfdp).' * lam / pert; %PfdpBeqz To 
    end
    %PfdpBeqv num_G74
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertBeqv, St_PertBeqv] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertBeqv evaluated in x+pert
        [dSf_dPfdp_PertBeqv, dSt_dPfdp_PertBeqv] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertBeqv
        %dAbr_dPfdpPertBeqv evaluated in x+pert        
        [dAf_dPfdp_PertBeqv, dAt_dPfdp_PertBeqv] = dAbr_dsh(dSf_dPfdp_PertBeqv, dSt_dPfdp_PertBeqv, Sf_PertBeqv, St_PertBeqv); %dAbr_dPfdpPertBeqv              
        %2nd Derivatives of Abr w.r.t. PfdpBeqv
        d2Af_dBeqvPfdp(:, k) = (dAf_dPfdp_PertBeqv - dAf_dPfdp).' * lam / pert; %PfdpBeqv From
        d2At_dBeqvPfdp(:, k) = (dAt_dPfdp_PertBeqv - dAt_dPfdp).' * lam / pert; %PfdpBeqv To 
    end
    %PfdpPfsh num_G85
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfsh, St_PertPfsh] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertPfsh evaluated in x+pert
        [dSf_dPfdp_PertPfsh, dSt_dPfdp_PertPfsh] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertPfsh
        %dAbr_dPfdpPertPfsh evaluated in x+pert        
        [dAf_dPfdp_PertPfsh, dAt_dPfdp_PertPfsh] = dAbr_dsh(dSf_dPfdp_PertPfsh, dSt_dPfdp_PertPfsh, Sf_PertPfsh, St_PertPfsh); %dAbr_dPfdpPertPfsh              
        %2nd Derivatives of Abr w.r.t. PfdpPfsh
        d2Af_dPfshPfdp(:, k) = (dAf_dPfdp_PertPfsh - dAf_dPfdp).' * lam / pert; %PfdpPfsh From
        d2At_dPfshPfdp(:, k) = (dAt_dPfdp_PertPfsh - dAt_dPfdp).' * lam / pert; %PfdpPfsh To 
    end    
    %PfdpQtma num_G86
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertQtma, St_PertQtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertQtma evaluated in x+pert
        [dSf_dPfdp_PertQtma, dSt_dPfdp_PertQtma] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertQtma
        %dAbr_dPfdpPertQtma evaluated in x+pert        
        [dAf_dPfdp_PertQtma, dAt_dPfdp_PertQtma] = dAbr_dsh(dSf_dPfdp_PertQtma, dSt_dPfdp_PertQtma, Sf_PertQtma, St_PertQtma); %dAbr_dPfdpPertPfsh              
        %2nd Derivatives of Abr w.r.t. PfdpQtma
        d2Af_dQtmaPfdp(:, k) = (dAf_dPfdp_PertQtma - dAf_dPfdp).' * lam / pert; %PfdpQtma From
        d2At_dQtmaPfdp(:, k) = (dAt_dPfdp_PertQtma - dAt_dPfdp).' * lam / pert; %PfdpQtma To 
    end 
    %PfdpVtma num_G87
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertVtma evaluated in x+pert
        [dSf_dPfdp_PertPfdp, dSt_dPfdp_PertPfdp] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertVtma
        %dAbr_dPfdpPertVtma evaluated in x+pert        
        [dAf_dPfdp_PertPfdp, dAt_dPfdp_PertPfdp] = dAbr_dsh(dSf_dPfdp_PertPfdp, dSt_dPfdp_PertPfdp, Sf_PertPfdp, St_PertPfdp); %dAbr_dPfdpPertPfsh              
        %2nd Derivatives of Abr w.r.t. PfdpVtma
        d2Af_dVtmaPfdp(:, k) = (dAf_dPfdp_PertPfdp - dAf_dPfdp).' * lam / pert; %PfdpVtma From
        d2At_dVtmaPfdp(:, k) = (dAt_dPfdp_PertPfdp - dAt_dPfdp).' * lam / pert; %PfdpVtma To 
    end 
    
    %VxVtma num_G18 num_G28
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertPfdp evaluated in x+pert        
        [dSf_dV1_PertPfdp, dSf_dV2_PertPfdp, dSt_dV1_PertPfdp, dSt_dV2_PertPfdp, Sf_PertPfdp, St_PertPfdp] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVxPertPfdp evaluated in x+pert
        [dAf_dV1_PertPfdp, dAf_dV2_PertPfdp, dAt_dV1_PertPfdp, dAt_dV2_PertPfdp] = dAbr_dV(dSf_dV1_PertPfdp, dSf_dV2_PertPfdp, dSt_dV1_PertPfdp, dSt_dV2_PertPfdp, Sf_PertPfdp, St_PertPfdp);
        %2nd Derivatives of Abr w.r.t. PfdpVx
        d2Af_dPfdpVa(:, k) = (dAf_dV1_PertPfdp - dAf_dV1).' * lam / pert;  %VaPfdp from, size of [nb, nPfdp] 
        d2Af_dPfdpVm(:, k) = (dAf_dV2_PertPfdp - dAf_dV2).' * lam / pert;  %VmPfdp from, size of [nb, nPfdp]
        d2At_dPfdpVa(:, k) = (dAt_dV1_PertPfdp - dAt_dV1).' * lam / pert;  %VaPfdp  to , size of [nb, nPfdp] 
        d2At_dPfdpVm(:, k) = (dAt_dV2_PertPfdp - dAt_dV2).' * lam / pert;  %VmPfdp  to , size of [nb, nPfdp]        
    end
    %BeqzPfdp num_G38
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dBeqzPertPfdp evaluated in x+pert
        [dSf_dBeqz_PertPfdp, dSt_dBeqz_PertPfdp] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertPfdp
        %dAbr_dBeqzPertPfdp evaluated in x+pert        
        [dAf_dBeqz_PertPfdp, dAt_dBeqz_PertPfdp] = dAbr_dBeq(dSf_dBeqz_PertPfdp, dSt_dBeqz_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Sbr w.r.t. BeqzPfdp
        d2Af_dPfdpBeqz(:, k) = (dAf_dBeqz_PertPfdp - dAf_dBeqz).' * lam / pert;  %BeqzPfdp from, size of [nBeqz, nPfdp]
        d2At_dPfdpBeqz(:, k) = (dAt_dBeqz_PertPfdp - dAt_dBeqz).' * lam / pert;  %BeqzPfdp  to , size of [nBeqz, nPfdp]        
     end
    %BeqvPfdp num_G48
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dBeqvPertPfdp evaluated in x+pert
        [dSf_dBeqv_PertPfdp, dSt_dBeqv_PertPfdp] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertPfdp
        %dAbr_dBeqvPertPfdp evaluated in x+pert        
        [dAf_dBeqv_PertPfdp, dAt_dBeqv_PertPfdp] = dAbr_dBeq(dSf_dBeqv_PertPfdp, dSt_dBeqv_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Sbr w.r.t. BeqvPfdp
        d2Af_dPfdpBeqv(:, k) = (dAf_dBeqv_PertPfdp - dAf_dBeqv).' * lam / pert;  %BeqvPfdp from, size of [nBeqv, nPfdp]
        d2At_dPfdpBeqv(:, k) = (dAt_dBeqv_PertPfdp - dAt_dBeqv).' * lam / pert;  %BeqvPfdp  to , size of [nBeqv, nPfdp]        
    end
    %PfshPfdp num_G58
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfshPertPfdp evaluated in x+pert
        [dSf_dPfsh_PertPfdp, dSt_dPfsh_PertPfdp] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertPfdp
        %dAbr_dPfshPertPfdp evaluated in x+pert        
        [dAf_dPfsh_PertPfdp, dAt_dPfsh_PertPfdp] = dAbr_dsh(dSf_dPfsh_PertPfdp, dSt_dPfsh_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Sbr w.r.t. PfshPfdp
        d2Af_dPfdpPfsh(:, k) = (dAf_dPfsh_PertPfdp - dAf_dPfsh).' * lam / pert;  %PfshPfdp from, size of [nPfsh, nPfdp]
        d2At_dPfdpPfsh(:, k) = (dAt_dPfsh_PertPfdp - dAt_dPfsh).' * lam / pert;  %PfshPfdp  to , size of [nPfsh, nPfdp]        
    end
    %QtmaPfdp num_G68
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dQtmaPertPfdp evaluated in x+pert
        [dSf_dQtma_PertPfdp, dSt_dQtma_PertPfdp] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertPfdp
        %dAbr_dQtmaPertPfdp evaluated in x+pert        
        [dAf_dQtma_PertPfdp, dAt_dQtma_PertPfdp] = dAbr_dma(dSf_dQtma_PertPfdp, dSt_dQtma_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Sbr w.r.t. QtmaPfdp
        d2Af_dPfdpQtma(:, k) = (dAf_dQtma_PertPfdp - dAf_dQtma).' * lam / pert;  %QtmaPfdp from, size of [nQtma, nPfdp]
        d2At_dPfdpQtma(:, k) = (dAt_dQtma_PertPfdp - dAt_dQtma).' * lam / pert;  %QtmaPfdp  to , size of [nQtma, nPfdp]        
    end
    %VtmaPfdp num_G78
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dVtmaPertPfdp evaluated in x+pert
        [dSf_dVtma_PertPfdp, dSt_dVtma_PertPfdp] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertPfdp
        %dAbr_dVtmaPertPfdp evaluated in x+pert        
        [dAf_dVtma_PertPfdp, dAt_dVtma_PertPfdp] = dAbr_dma(dSf_dVtma_PertPfdp, dSt_dVtma_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Sbr w.r.t. VtmaPfdp
        d2Af_dPfdpVtma(:, k) = (dAf_dVtma_PertPfdp - dAf_dVtma).' * lam / pert;  %VtmaPfdp from, size of [nVtma, nPfdp]
        d2At_dPfdpVtma(:, k) = (dAt_dVtma_PertPfdp - dAt_dVtma).' * lam / pert;  %VtmaPfdp  to , size of [nVtma, nPfdp]        
    end    
    %Pfdp2 num_G88
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagPfdpsel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfdp, St_PertPfdp] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfdpPertPfdp evaluated in x+pert
        [dSf_dPfdp_PertPfdp, dSt_dPfdp_PertPfdp] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertPfdp
        %dAbr_dPfdpPertPfdp evaluated in x+pert        
        [dAf_dPfdp_PertPfdp, dAt_dPfdp_PertPfdp] = dAbr_dsh(dSf_dPfdp_PertPfdp, dSt_dPfdp_PertPfdp, Sf_PertPfdp, St_PertPfdp);             
        %2nd Derivatives of Abr w.r.t. Pfdp2   
        d2Af_dPfdp2(:, k) = (dAf_dPfdp_PertPfdp - dAf_dPfdp).' * lam / pert;  %Pfdp2 from, size of [nPfdp , nPfdp] 
        d2At_dPfdp2(:, k) = (dAt_dPfdp_PertPfdp - dAt_dPfdp).' * lam / pert;  %Pfdp2  to , size of [nPfdp , nPfdp]  
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf18 = sparse(d2Af_dPfdpVa);
num_Hf28 = sparse(d2Af_dPfdpVm);
num_Hf38 = sparse(d2Af_dPfdpBeqz);
num_Hf48 = sparse(d2Af_dPfdpBeqv);
num_Hf58 = sparse(d2Af_dPfdpPfsh);
num_Hf68 = sparse(d2Af_dPfdpQtma);
num_Hf78 = sparse(d2Af_dPfdpVtma);

num_Hf81 = sparse(d2Af_dVaPfdp);
num_Hf82 = sparse(d2Af_dVmPfdp);
num_Hf83 = sparse(d2Af_dBeqzPfdp);
num_Hf84 = sparse(d2Af_dBeqvPfdp);
num_Hf85 = sparse(d2Af_dPfshPfdp);
num_Hf86 = sparse(d2Af_dQtmaPfdp);
num_Hf87 = sparse(d2Af_dVtmaPfdp);

num_Hf88 = sparse(d2Af_dPfdp2);

num_Ht18 = sparse(d2At_dPfdpVa);
num_Ht28 = sparse(d2At_dPfdpVm);
num_Ht38 = sparse(d2At_dPfdpBeqz);
num_Ht48 = sparse(d2At_dPfdpBeqv);
num_Ht58 = sparse(d2At_dPfdpPfsh);
num_Ht68 = sparse(d2At_dPfdpQtma);
num_Ht78 = sparse(d2At_dPfdpVtma);

num_Ht81 = sparse(d2At_dVaPfdp);
num_Ht82 = sparse(d2At_dVmPfdp);
num_Ht83 = sparse(d2At_dBeqzPfdp);
num_Ht84 = sparse(d2At_dBeqvPfdp);
num_Ht85 = sparse(d2At_dPfshPfdp);
num_Ht86 = sparse(d2At_dQtmaPfdp);
num_Ht87 = sparse(d2At_dVtmaPfdp);

num_Ht88 = sparse(d2At_dPfdp2);
