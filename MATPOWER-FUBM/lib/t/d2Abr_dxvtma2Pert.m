function [num_Hf17, num_Hf27, num_Hf37, num_Hf47, num_Hf57, num_Hf67, num_Hf71, num_Hf72, num_Hf73, num_Hf74, num_Hf75, num_Hf76, num_Hf77,...
          num_Ht17, num_Ht27, num_Ht37, num_Ht47, num_Ht57, num_Ht67, num_Ht71, num_Ht72, num_Ht73, num_Ht74, num_Ht75, num_Ht76, num_Ht77] = d2Abr_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXVTMA2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma, vtmavtma(Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%  [NUM_HF17, NUM_HF27, NUM_HF37, NUM_HF47, NUM_HF57, NUM_HF67, NUM_HF71, NUM_HF72, NUM_HF73, NUM_HF74, NUM_HF75, NUM_HF76, NUM_HF77,...
%   NUM_HT17, NUM_HT27, NUM_HT37, NUM_HT47, NUM_HT57, NUM_HT67, NUM_HT71, NUM_HT72, NUM_HT73, NUM_HT74, NUM_HT75, NUM_HT76, NUM_HT77] = D2ABR_DXVTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 26 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H17 = HvtmaVa   = d2Abr_dvtmaVa     = d/dvtma ( (dAbr/dVa    )'*mu )
%   H27 = HvtmaVm   = d2Abr_dvtmaVm     = d/dvtma ( (dAbr/dVm    )'*mu )
%   H37 = HvtmaBz   = d2Abr_dvtmaBeqz   = d/dvtma ( (dAbr/dBeqz  )'*mu )
%   H47 = HvtmaBv   = d2Abr_dvtmaBeqv   = d/dvtma ( (dAbr/dBeqv  )'*mu )
%   H57 = HvtmaSh   = d2Abr_dvtmaSh     = d/dvtma ( (dAbr/dSh    )'*mu )
%   H57 = HvtmaSh   = d2Abr_dvtmaqtma   = d/dvtma ( (dAbr/dqtma  )'*mu )
%   H71 = HVavtma   = d2Abr_dVavtma     = d/dVa   ( (dAbr/dvtma  )'*mu )
%   H72 = HVmvtma   = d2Abr_dVmvtma     = d/dVm   ( (dAbr/dvtma  )'*mu )
%   H73 = HBzvtma   = d2Abr_dBeqzvtma   = d/dBeqz ( (dAbr/dvtma  )'*mu )
%   H74 = HBvvtma   = d2Abr_dBeqvvtma   = d/dBeqv ( (dAbr/dvtma  )'*mu )
%   H75 = HShvtma   = d2Abr_dShvtma     = d/dsh   ( (dAbr/dvtma  )'*mu )
%   H76 = Hqtmavtma = d2Abr_dqtmavtma   = d/dqtma ( (dAbr/dvtma  )'*mu )
%   H77 = Hvtmavtma = d2Abr_dvtmavtma   = d/dvtma ( (dAbr/dvtma  )'*mu )
%   
%   H17 = HvtmaVa   = d2Abr_dvtmaVa     =  (  2*real(d2Sbr_dvtmaVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dvtma   ))'*mu 
%   H27 = HvtmaVm   = d2Abr_dvtmaVm     =  (  2*real(d2Sbr_dvtmaVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dvtma   ))'*mu
%   H37 = HvtmaBz   = d2Abr_dvtmaBeqz   =  (  2*real(d2Sbr_dvtmaBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dvtma   ))'*mu
%   H47 = HvtmaBv   = d2Abr_dvtmaBeqv   =  (  2*real(d2Sbr_dvtmaBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dvtma   ))'*mu
%   H57 = HvtmaSh   = d2Abr_dvtmaSh     =  (  2*real(d2Sbr_dvtmaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dvtma   ))'*mu
%   H67 = Hvtmaqtma = d2Abr_dvtmaqtma   =  (  2*real(d2Sbr_dvtmaqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dvtma   ))'*mu
%   H71 = HVavtma   = d2Abr_dVavtma     =  (  2*real(d2Sbr_dVavtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dVa     ))'*mu
%   H72 = HVmvtma   = d2Abr_dVmvtma     =  (  2*real(d2Sbr_dVmvtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dVm     ))'*mu
%   H73 = HBzvtma   = d2Abr_dBeqzvtma   =  (  2*real(d2Sbr_dBeqzvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dBeqz   ))'*mu
%   H74 = HBvvtma   = d2Abr_dBeqvvtma   =  (  2*real(d2Sbr_dBeqvvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dBeqv   ))'*mu
%   H75 = HShvtma   = d2Abr_dShvtma     =  (  2*real(d2Sbr_dShvtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dSh     ))'*mu 
%   H76 = Hqtmavtma = d2Abr_dqtmavtma   =  (  2*real(d2Sbr_dqtmavtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dqtma   ))'*mu 
%   H77 = Hvtmavtma = d2Abr_dvtmavtma   =  (  2*real(d2Sbr_dvtmavtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dvtma   ))'*mu 
%
%   Example:
%  [NUM_HF17, NUM_HF27, NUM_HF37, NUM_HF47, NUM_HF57, NUM_HF67, NUM_HF71, NUM_HF72, NUM_HF73, NUM_HF74, NUM_HF75, NUM_HF76, NUM_HF77,...
%   NUM_HT17, NUM_HT27, NUM_HT37, NUM_HT47, NUM_HT57, NUM_HT67, NUM_HT71, NUM_HT72, NUM_HT73, NUM_HT74, NUM_HT75, NUM_HT76, NUM_HT77] = D2ABR_DXVTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
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
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

%% Calculation of derivatives
if vcart
    error('d2Abr_dxvtma2Pert: Derivatives of Flow Limit equations w.r.t ma using Finite Differences Method in cartasian have not been coded yet')    

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
    
    %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St); 
    [dAf_dPfsh, dAt_dPfsh] = dAbr_dsh(dSf_dPfsh, dSt_dPfsh, Sf, St); 
    [dAf_dQtma, dAt_dQtma] = dAbr_dma(dSf_dQtma, dSt_dQtma, Sf, St); 
    [dAf_dVtma, dAt_dVtma] = dAbr_dma(dSf_dVtma, dSt_dVtma, Sf, St); 
    
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
     
    %Dimensionalize (Allocate for computational speed)  
    d2Af_dVtmaVa   = sparse( zeros(nb,   nVtma) ); 
    d2Af_dVtmaVm   = sparse( zeros(nb,   nVtma) );
    d2Af_dVtmaBeqz = sparse( zeros(nBeqz,nVtma) );
    d2Af_dVtmaBeqv = sparse( zeros(nBeqv,nVtma) );
    d2Af_dVtmaPfsh = sparse( zeros(nPfsh,nVtma) );   
    d2Af_dVtmaQtma = sparse( zeros(nQtma,nVtma) );  
    
    d2Af_dVaVtma   = sparse( zeros(nVtma,nb   ) ); 
    d2Af_dVmVtma   = sparse( zeros(nVtma,nb   ) );
    d2Af_dBeqzVtma = sparse( zeros(nVtma,nBeqz) );
    d2Af_dBeqvVtma = sparse( zeros(nVtma,nBeqv) );
    d2Af_dPfshVtma = sparse( zeros(nVtma,nPfsh) );   
    d2Af_dQtmaVtma = sparse( zeros(nVtma,nQtma) );     
    
    d2Af_dVtma2    = sparse( zeros(nVtma,nVtma) );

    d2At_dVtmaVa   = sparse( zeros(nb,   nVtma) ); 
    d2At_dVtmaVm   = sparse( zeros(nb,   nVtma) );
    d2At_dVtmaBeqz = sparse( zeros(nBeqz,nVtma) );
    d2At_dVtmaBeqv = sparse( zeros(nBeqv,nVtma) );
    d2At_dVtmaPfsh = sparse( zeros(nPfsh,nVtma) );   
    d2At_dVtmaQtma = sparse( zeros(nQtma,nVtma) );  
    
    d2At_dVaVtma   = sparse( zeros(nVtma,nb   ) ); 
    d2At_dVmVtma   = sparse( zeros(nVtma,nb   ) );
    d2At_dBeqzVtma = sparse( zeros(nVtma,nBeqz) );
    d2At_dBeqvVtma = sparse( zeros(nVtma,nBeqv) );
    d2At_dPfshVtma = sparse( zeros(nVtma,nPfsh) );   
    d2At_dQtmaVtma = sparse( zeros(nVtma,nQtma) );     
   
    d2At_dVtma2    = sparse( zeros(nVtma,nVtma) );
    
    %VtmaVa num_G71
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dVtma_PertVa, dSt_dVtma_PertVa] = dSbr_dma(branch, V1p, 4, vcart); %dSbr_dVtmaPertVa
        [dAf_dVtma_PertVa, dAt_dVtma_PertVa] = dAbr_dma(dSf_dVtma_PertVa, dSt_dVtma_PertVa, Sf_PertVa, St_PertVa); %dAbr_dVtmaPertVa       
        d2Af_dVaVtma(:, k) = (dAf_dVtma_PertVa - dAf_dVtma).' * lam / pert; %VtmaVa from, size of [nVtma, nb]
        d2At_dVaVtma(:, k) = (dAt_dVtma_PertVa - dAt_dVtma).' * lam / pert; %VtmaVa  to , size of [nVtma, nb]
    end
    %VtmaVm num_G72
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dVtma_PertVm, dSt_dVtma_PertVm] = dSbr_dma(branch, V2p, 4, vcart); %dSbr_dVtmaPertVm
        [dAf_dVtma_PertVm, dAt_dVtma_PertVm] = dAbr_dma(dSf_dVtma_PertVm, dSt_dVtma_PertVm, Sf_PertVm, St_PertVm); %dAbr_dVtmaPertVm       
        d2Af_dVmVtma(:, k) = (dAf_dVtma_PertVm - dAf_dVtma).' * lam / pert; %VtmaVm from, size of [nVtma, nb]
        d2At_dVmVtma(:, k) = (dAt_dVtma_PertVm - dAt_dVtma).' * lam / pert; %VtmaVm  to , size of [nVtma, nb]
    end
    %VtmaBeqz num_G73
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
        %dSbr_dVtmaPertBeqz evaluated in x+pert
        [dSf_dVtma_PertBeqz, dSt_dVtma_PertBeqz] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertBeqz
        %dAbr_dVtmaPertBeqz evaluated in x+pert        
        [dAf_dVtma_PertBeqz, dAt_dVtma_PertBeqz] = dAbr_dma(dSf_dVtma_PertBeqz, dSt_dVtma_PertBeqz, Sf_PertBeqz, St_PertBeqz); %dAbr_dVtmaPertBeqz              
        %2nd Derivatives of Abr w.r.t. VtmaBeqz
        d2Af_dBeqzVtma(:, k) = (dAf_dVtma_PertBeqz - dAf_dVtma).' * lam / pert; %VtmaBeqz From
        d2At_dBeqzVtma(:, k) = (dAt_dVtma_PertBeqz - dAt_dVtma).' * lam / pert; %VtmaBeqz To 
    end
    %VtmaBeqv num_G74
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
        %dSbr_dVtmaPertBeqv evaluated in x+pert
        [dSf_dVtma_PertBeqv, dSt_dVtma_PertBeqv] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertBeqv
        %dAbr_dVtmaPertBeqv evaluated in x+pert        
        [dAf_dVtma_PertBeqv, dAt_dVtma_PertBeqv] = dAbr_dma(dSf_dVtma_PertBeqv, dSt_dVtma_PertBeqv, Sf_PertBeqv, St_PertBeqv); %dAbr_dVtmaPertBeqv              
        %2nd Derivatives of Abr w.r.t. VtmaBeqv
        d2Af_dBeqvVtma(:, k) = (dAf_dVtma_PertBeqv - dAf_dVtma).' * lam / pert; %VtmaBeqv From
        d2At_dBeqvVtma(:, k) = (dAt_dVtma_PertBeqv - dAt_dVtma).' * lam / pert; %VtmaBeqv To 
    end
    %VtmaPfsh num_G75
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
        %dSbr_dVtmaPertPfsh evaluated in x+pert
        [dSf_dVtma_PertPfsh, dSt_dVtma_PertPfsh] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertPfsh
        %dAbr_dVtmaPertPfsh evaluated in x+pert        
        [dAf_dVtma_PertPfsh, dAt_dVtma_PertPfsh] = dAbr_dma(dSf_dVtma_PertPfsh, dSt_dVtma_PertPfsh, Sf_PertPfsh, St_PertPfsh); %dAbr_dVtmaPertPfsh              
        %2nd Derivatives of Abr w.r.t. VtmaPfsh
        d2Af_dPfshVtma(:, k) = (dAf_dVtma_PertPfsh - dAf_dVtma).' * lam / pert; %VtmaPfsh From
        d2At_dPfshVtma(:, k) = (dAt_dVtma_PertPfsh - dAt_dVtma).' * lam / pert; %VtmaPfsh To 
    end    
    %VtmaQtma num_G76
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
        %dSbr_dVtmaPertQtma evaluated in x+pert
        [dSf_dVtma_PertQtma, dSt_dVtma_PertQtma] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertQtma
        %dAbr_dVtmaPertQtma evaluated in x+pert        
        [dAf_dVtma_PertQtma, dAt_dVtma_PertQtma] = dAbr_dma(dSf_dVtma_PertQtma, dSt_dVtma_PertQtma, Sf_PertQtma, St_PertQtma); %dAbr_dVtmaPertPfsh              
        %2nd Derivatives of Abr w.r.t. VtmaQtma
        d2Af_dQtmaVtma(:, k) = (dAf_dVtma_PertQtma - dAf_dVtma).' * lam / pert; %VtmaQtma From
        d2At_dQtmaVtma(:, k) = (dAt_dVtma_PertQtma - dAt_dVtma).' * lam / pert; %VtmaQtma To 
    end 
    
    %VxVtma num_G17 num_G27
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertVtma evaluated in x+pert        
        [dSf_dV1_PertVtma, dSf_dV2_PertVtma, dSt_dV1_PertVtma, dSt_dV2_PertVtma, Sf_PertVtma, St_PertVtma] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVxPertVtma evaluated in x+pert
        [dAf_dV1_PertVtma, dAf_dV2_PertVtma, dAt_dV1_PertVtma, dAt_dV2_PertVtma] = dAbr_dV(dSf_dV1_PertVtma, dSf_dV2_PertVtma, dSt_dV1_PertVtma, dSt_dV2_PertVtma, Sf_PertVtma, St_PertVtma);
        %2nd Derivatives of Abr w.r.t. VtmaVx
        d2Af_dVtmaVa(:, k) = (dAf_dV1_PertVtma - dAf_dV1).' * lam / pert;  %VaVtma from, size of [nb, nVtma] 
        d2Af_dVtmaVm(:, k) = (dAf_dV2_PertVtma - dAf_dV2).' * lam / pert;  %VmVtma from, size of [nb, nVtma]
        d2At_dVtmaVa(:, k) = (dAt_dV1_PertVtma - dAt_dV1).' * lam / pert;  %VaVtma  to , size of [nb, nVtma] 
        d2At_dVtmaVm(:, k) = (dAt_dV2_PertVtma - dAt_dV2).' * lam / pert;  %VmVtma  to , size of [nb, nVtma]        
 end
    %BeqzVtma num_G37
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertVtma, St_PertVtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dBeqzPertVtma evaluated in x+pert
        [dSf_dBeqz_PertVtma, dSt_dBeqz_PertVtma] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertVtma
        %dAbr_dBeqzPertVtma evaluated in x+pert        
        [dAf_dBeqz_PertVtma, dAt_dBeqz_PertVtma] = dAbr_dBeq(dSf_dBeqz_PertVtma, dSt_dBeqz_PertVtma, Sf_PertVtma, St_PertVtma);             
        %2nd Derivatives of Sbr w.r.t. BeqzVtma
        d2Af_dVtmaBeqz(:, k) = (dAf_dBeqz_PertVtma - dAf_dBeqz).' * lam / pert;  %BeqzVtma from, size of [nBeqz, nVtma]
        d2At_dVtmaBeqz(:, k) = (dAt_dBeqz_PertVtma - dAt_dBeqz).' * lam / pert;  %BeqzVtma  to , size of [nBeqz, nVtma]        
     end
    %BeqvVtma num_G47
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertVtma, St_PertVtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dBeqvPertVtma evaluated in x+pert
        [dSf_dBeqv_PertVtma, dSt_dBeqv_PertVtma] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertVtma
        %dAbr_dBeqvPertVtma evaluated in x+pert        
        [dAf_dBeqv_PertVtma, dAt_dBeqv_PertVtma] = dAbr_dBeq(dSf_dBeqv_PertVtma, dSt_dBeqv_PertVtma, Sf_PertVtma, St_PertVtma);             
        %2nd Derivatives of Sbr w.r.t. BeqvVtma
        d2Af_dVtmaBeqv(:, k) = (dAf_dBeqv_PertVtma - dAf_dBeqv).' * lam / pert;  %BeqvVtma from, size of [nBeqv, nVtma]
        d2At_dVtmaBeqv(:, k) = (dAt_dBeqv_PertVtma - dAt_dBeqv).' * lam / pert;  %BeqvVtma  to , size of [nBeqv, nVtma]        
    end
    %PfshVtma num_G57
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertVtma, St_PertVtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfshPertVtma evaluated in x+pert
        [dSf_dPfsh_PertVtma, dSt_dPfsh_PertVtma] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertVtma
        %dAbr_dPfshPertVtma evaluated in x+pert        
        [dAf_dPfsh_PertVtma, dAt_dPfsh_PertVtma] = dAbr_dsh(dSf_dPfsh_PertVtma, dSt_dPfsh_PertVtma, Sf_PertVtma, St_PertVtma);             
        %2nd Derivatives of Sbr w.r.t. PfshVtma
        d2Af_dVtmaPfsh(:, k) = (dAf_dPfsh_PertVtma - dAf_dPfsh).' * lam / pert;  %PfshVtma from, size of [nPfsh, nVtma]
        d2At_dVtmaPfsh(:, k) = (dAt_dPfsh_PertVtma - dAt_dPfsh).' * lam / pert;  %PfshVtma  to , size of [nPfsh, nVtma]        
    end
    %QtmaVtma num_G67
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertVtma, St_PertVtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dQtmaPertVtma evaluated in x+pert
        [dSf_dQtma_PertVtma, dSt_dQtma_PertVtma] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertVtma
        %dAbr_dQtmaPertVtma evaluated in x+pert        
        [dAf_dQtma_PertVtma, dAt_dQtma_PertVtma] = dAbr_dma(dSf_dQtma_PertVtma, dSt_dQtma_PertVtma, Sf_PertVtma, St_PertVtma);             
        %2nd Derivatives of Sbr w.r.t. QtmaVtma
        d2Af_dVtmaQtma(:, k) = (dAf_dQtma_PertVtma - dAf_dQtma).' * lam / pert;  %QtmaVtma from, size of [nQtma, nVtma]
        d2At_dVtmaQtma(:, k) = (dAt_dQtma_PertVtma - dAt_dQtma).' * lam / pert;  %QtmaVtma  to , size of [nQtma, nVtma]        
    end
    %Vtma2 num_G77
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertVtma, St_PertVtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dVtmaPertVtma evaluated in x+pert
        [dSf_dVtma_PertVtma, dSt_dVtma_PertVtma] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertVtma
        %dAbr_dVtmaPertVtma evaluated in x+pert        
        [dAf_dVtma_PertVtma, dAt_dVtma_PertVtma] = dAbr_dma(dSf_dVtma_PertVtma, dSt_dVtma_PertVtma, Sf_PertVtma, St_PertVtma);             
        %2nd Derivatives of Abr w.r.t. Vtma2   
        d2Af_dVtma2(:, k) = (dAf_dVtma_PertVtma - dAf_dVtma).' * lam / pert;  %Vtma2 from, size of [nVtma , nVtma] 
        d2At_dVtma2(:, k) = (dAt_dVtma_PertVtma - dAt_dVtma).' * lam / pert;  %Vtma2  to , size of [nVtma , nVtma]  
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf17 = sparse(d2Af_dVtmaVa);
num_Hf27 = sparse(d2Af_dVtmaVm);
num_Hf37 = sparse(d2Af_dVtmaBeqz);
num_Hf47 = sparse(d2Af_dVtmaBeqv);
num_Hf57 = sparse(d2Af_dVtmaPfsh);
num_Hf67 = sparse(d2Af_dVtmaQtma);

num_Hf71 = sparse(d2Af_dVaVtma);
num_Hf72 = sparse(d2Af_dVmVtma);
num_Hf73 = sparse(d2Af_dBeqzVtma);
num_Hf74 = sparse(d2Af_dBeqvVtma);
num_Hf75 = sparse(d2Af_dPfshVtma);
num_Hf76 = sparse(d2Af_dQtmaVtma);

num_Hf77 = sparse(d2Af_dVtma2);

num_Ht17 = sparse(d2At_dVtmaVa);
num_Ht27 = sparse(d2At_dVtmaVm);
num_Ht37 = sparse(d2At_dVtmaBeqz);
num_Ht47 = sparse(d2At_dVtmaBeqv);
num_Ht57 = sparse(d2At_dVtmaPfsh);
num_Ht67 = sparse(d2At_dVtmaQtma);

num_Ht71 = sparse(d2At_dVaVtma);
num_Ht72 = sparse(d2At_dVmVtma);
num_Ht73 = sparse(d2At_dBeqzVtma);
num_Ht74 = sparse(d2At_dBeqvVtma);
num_Ht75 = sparse(d2At_dPfshVtma);
num_Ht76 = sparse(d2At_dQtmaVtma);

num_Ht77 = sparse(d2At_dVtma2);
