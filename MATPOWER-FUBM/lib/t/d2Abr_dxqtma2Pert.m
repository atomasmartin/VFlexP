function [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66,...
          num_Ht16, num_Ht26, num_Ht36, num_Ht46, num_Ht56, num_Ht61, num_Ht62, num_Ht63, num_Ht64, num_Ht65, num_Ht66] = d2Abr_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXQTMA2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [NUM_HF16, NUM_HF26, NUM_HF36, NUM_HF46, NUM_HF56, NUM_HF61, NUM_HF62, NUM_HF63, NUM_HF64, NUM_HF65, NUM_HF66,...
%    NUM_HT16, NUM_HT26, NUM_HT36, NUM_HT46, NUM_HT56, NUM_HT61, NUM_HT62, NUM_HT63, NUM_HT64, NUM_HT65, NUM_HT66] = D2ABR_DXQTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 22 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H16 = HqtmaVa   = d2Abr_dqtmaVa     = d/dqtma ( (dAbr/dVa    )'*mu )
%   H26 = HqtmaVm   = d2Abr_dqtmaVm     = d/dqtma ( (dAbr/dVm    )'*mu )
%   H36 = HqtmaBz   = d2Abr_dqtmaBeqz   = d/dqtma ( (dAbr/dBeqz  )'*mu )
%   H46 = HqtmaBv   = d2Abr_dqtmaBeqv   = d/dqtma ( (dAbr/dBeqv  )'*mu )
%   H56 = HqtmaSh   = d2Abr_dqtmaSh     = d/dqtma ( (dAbr/dSh    )'*mu )
%   H61 = HVaqtma   = d2Abr_dVaqtma     = d/dVa   ( (dAbr/dqtma  )'*mu )
%   H62 = HVmqtma   = d2Abr_dVmqtma     = d/dVm   ( (dAbr/dqtma  )'*mu )
%   H63 = HBzqtma   = d2Abr_dBeqzqtma   = d/dBeqz ( (dAbr/dqtma  )'*mu )
%   H64 = HBvqtma   = d2Abr_dBeqvqtma   = d/dBeqv ( (dAbr/dqtma  )'*mu )
%   H65 = HShqtma   = d2Abr_dShqtma     = d/dsh   ( (dAbr/dqtma  )'*mu )
%   H66 = Hqtmaqtma = d2Abr_dqtmaqtma   = d/dqtma ( (dAbr/dqtma  )'*mu )
%   
%   H16 = HqtmaVa   = d2Abr_dqtmaVa     =  (  2*real(d2Sbr_dqtmaVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dqtma   ))'*mu 
%   H26 = HqtmaVm   = d2Abr_dqtmaVm     =  (  2*real(d2Sbr_dqtmaVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dqtma   ))'*mu
%   H36 = HqtmaBz   = d2Abr_dqtmaBeqz   =  (  2*real(d2Sbr_dqtmaBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dqtma   ))'*mu
%   H46 = HqtmaBv   = d2Abr_dqtmaBeqv   =  (  2*real(d2Sbr_dqtmaBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dqtma   ))'*mu
%   H56 = HqtmaSh   = d2Abr_dqtmaSh     =  (  2*real(d2Sbr_dqtmaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dqtma   ))'*mu
%   H61 = HVaqtma   = d2Abr_dVaqtma     =  (  2*real(d2Sbr_dVaqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dVa     ))'*mu
%   H62 = HVmqtma   = d2Abr_dVmqtma     =  (  2*real(d2Sbr_dVmqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dVm     ))'*mu
%   H63 = HBzqtma   = d2Abr_dBeqzqtma   =  (  2*real(d2Sbr_dBeqzqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dBeqz   ))'*mu
%   H64 = HBvqtma   = d2Abr_dBeqvqtma   =  (  2*real(d2Sbr_dBeqvqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dBeqv   ))'*mu
%   H65 = HShqtma   = d2Abr_dShqtma     =  (  2*real(d2Sbr_dShqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dSh     ))'*mu 
%   H66 = Hqtmaqtma = d2Abr_dqtmaqtma   =  (  2*real(d2Sbr_dqtmaqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dma     ))'*mu 
%
%   Example:
%   [NUM_HF16, NUM_HF26, NUM_HF36, NUM_HF46, NUM_HF56, NUM_HF61, NUM_HF62, NUM_HF63, NUM_HF64, NUM_HF65, NUM_HF66,...
%    NUM_HT16, NUM_HT26, NUM_HT36, NUM_HT46, NUM_HT56, NUM_HT61, NUM_HT62, NUM_HT63, NUM_HT64, NUM_HT65, NUM_HT66] = D2ABR_DXQTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
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

%% Calculation of derivatives
if vcart
    error('d2Abr_dqtma2Pert: Derivatives of Flow Limit equations w.r.t ma using Finite Differences Method in cartasian have not been coded yet')    

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

    %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St); 
    [dAf_dPfsh, dAt_dPfsh] = dAbr_dsh(dSf_dPfsh, dSt_dPfsh, Sf, St); 
    [dAf_dQtma, dAt_dQtma] = dAbr_dma(dSf_dQtma, dSt_dQtma, Sf, St); 
 
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
    
    %Dimensionalize (Allocate for computational speed)  
    d2Af_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Af_dQtmaVm   = sparse( zeros(nb,   nQtma) );
    d2Af_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Af_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Af_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) ); 

    d2Af_dVaQtma   = sparse( zeros(nQtma,nb   ) ); 
    d2Af_dVmQtma   = sparse( zeros(nQtma,nb   ) );
    d2Af_dBeqzQtma = sparse( zeros(nQtma,nBeqz) );
    d2Af_dBeqvQtma = sparse( zeros(nQtma,nBeqv) );
    d2Af_dPfshQtma = sparse( zeros(nQtma,nPfsh) ); 
    
    d2Af_dQtma2    = sparse( zeros(nQtma,nQtma) );

    d2At_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2At_dQtmaVm   = sparse( zeros(nb,   nQtma) );
    d2At_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2At_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2At_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) ); 

    d2At_dVaQtma   = sparse( zeros(nQtma,nb   ) ); 
    d2At_dVmQtma   = sparse( zeros(nQtma,nb   ) );
    d2At_dBeqzQtma = sparse( zeros(nQtma,nBeqz) );
    d2At_dBeqvQtma = sparse( zeros(nQtma,nBeqv) );
    d2At_dPfshQtma = sparse( zeros(nQtma,nPfsh) ); 
    
    d2At_dQtma2    = sparse( zeros(nQtma,nQtma) );    
    
    %QtmaVa num_G61
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dQtma_PertVa, dSt_dQtma_PertVa] = dSbr_dma(branch, V1p, 2, vcart); %dSbr_dQtmaPertVa
        [dAf_dQtma_PertVa, dAt_dQtma_PertVa] = dAbr_dma(dSf_dQtma_PertVa, dSt_dQtma_PertVa, Sf_PertVa, St_PertVa); %dAbr_dQtmaPertVa       
        d2Af_dVaQtma(:, k) = (dAf_dQtma_PertVa - dAf_dQtma).' * lam / pert; %QtmaVa from, size of [nQtma, nb]
        d2At_dVaQtma(:, k) = (dAt_dQtma_PertVa - dAt_dQtma).' * lam / pert; %QtmaVa  to , size of [nQtma, nb]
    end
    %QtmaVm num_G62
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dQtma_PertVm, dSt_dQtma_PertVm] = dSbr_dma(branch, V2p, 2, vcart); %dSbr_dQtmaPertVm
        [dAf_dQtma_PertVm, dAt_dQtma_PertVm] = dAbr_dma(dSf_dQtma_PertVm, dSt_dQtma_PertVm, Sf_PertVm, St_PertVm); %dAbr_dQtmaPertVm       
        d2Af_dVmQtma(:, k) = (dAf_dQtma_PertVm - dAf_dQtma).' * lam / pert; %QtmaVm from, size of [nQtma, nb]
        d2At_dVmQtma(:, k) = (dAt_dQtma_PertVm - dAt_dQtma).' * lam / pert; %QtmaVm  to , size of [nQtma, nb]
    end
    %QtmaBeqz num_G63
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
        %dSbr_dQtmaPertBeqz evaluated in x+pert
        [dSf_dQtma_PertBeqz, dSt_dQtma_PertBeqz] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertBeqz
        %dAbr_dQtmaPertBeqz evaluated in x+pert        
        [dAf_dQtma_PertBeqz, dAt_dQtma_PertBeqz] = dAbr_dma(dSf_dQtma_PertBeqz, dSt_dQtma_PertBeqz, Sf_PertBeqz, St_PertBeqz); %dAbr_dQtmaPertBeqz              
        %2nd Derivatives of Abr w.r.t. QtmaBeqz
        d2Af_dBeqzQtma(:, k) = (dAf_dQtma_PertBeqz - dAf_dQtma).' * lam / pert; %QtmaBeqz From
        d2At_dBeqzQtma(:, k) = (dAt_dQtma_PertBeqz - dAt_dQtma).' * lam / pert; %QtmaBeqz To 
    end
    %QtmaBeqv num_G64
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
        %dSbr_dQtmaPertBeqv evaluated in x+pert
        [dSf_dQtma_PertBeqv, dSt_dQtma_PertBeqv] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertBeqv
        %dAbr_dPfshPertBeqv evaluated in x+pert        
        [dAf_dQtma_PertBeqv, dAt_dQtma_PertBeqv] = dAbr_dma(dSf_dQtma_PertBeqv, dSt_dQtma_PertBeqv, Sf_PertBeqv, St_PertBeqv); %dAbr_dQtmaPertBeqv              
        %2nd Derivatives of Abr w.r.t. QtmaBeqv
        d2Af_dBeqvQtma(:, k) = (dAf_dQtma_PertBeqv - dAf_dQtma).' * lam / pert; %QtmaBeqv From
        d2At_dBeqvQtma(:, k) = (dAt_dQtma_PertBeqv - dAt_dQtma).' * lam / pert; %QtmaBeqv To 
    end
    %QtmaPfsh num_G65
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
        %dSbr_dQtmaPertPfsh evaluated in x+pert
        [dSf_dQtma_PertPfsh, dSt_dQtma_PertPfsh] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertPfsh
        %dAbr_dPfshPertPfsh evaluated in x+pert        
        [dAf_dQtma_PertPfsh, dAt_dQtma_PertPfsh] = dAbr_dma(dSf_dQtma_PertPfsh, dSt_dQtma_PertPfsh, Sf_PertPfsh, St_PertPfsh); %dAbr_dQtmaPertPfsh              
        %2nd Derivatives of Abr w.r.t. QtmaPfsh
        d2Af_dPfshQtma(:, k) = (dAf_dQtma_PertPfsh - dAf_dQtma).' * lam / pert; %QtmaPfsh From
        d2At_dPfshQtma(:, k) = (dAt_dQtma_PertPfsh - dAt_dQtma).' * lam / pert; %QtmaPfsh To 
    end    
    
    %VxQtma num_G16 num_G26
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertQtma evaluated in x+pert        
        [dSf_dV1_PertQtma, dSf_dV2_PertQtma, dSt_dV1_PertQtma, dSt_dV2_PertQtma, Sf_PertQtma, St_PertQtma] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVxPertQtma evaluated in x+pert
        [dAf_dV1_PertQtma, dAf_dV2_PertQtma, dAt_dV1_PertQtma, dAt_dV2_PertQtma] = dAbr_dV(dSf_dV1_PertQtma, dSf_dV2_PertQtma, dSt_dV1_PertQtma, dSt_dV2_PertQtma, Sf_PertQtma, St_PertQtma);
        %2nd Derivatives of Abr w.r.t. QtmaVx
        d2Af_dQtmaVa(:, k) = (dAf_dV1_PertQtma - dAf_dV1).' * lam / pert;  %VaQtma from, size of [nb, nQtma] 
        d2Af_dQtmaVm(:, k) = (dAf_dV2_PertQtma - dAf_dV2).' * lam / pert;  %VmQtma from, size of [nb, nQtma]
        d2At_dQtmaVa(:, k) = (dAt_dV1_PertQtma - dAt_dV1).' * lam / pert;  %VaQtma  to , size of [nb, nQtma] 
        d2At_dQtmaVm(:, k) = (dAt_dV2_PertQtma - dAt_dV2).' * lam / pert;  %VmQtma  to , size of [nb, nQtma]        
 end
    %BeqzQtma num_G36
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
        %dSbr_dBeqzPertQtma evaluated in x+pert
        [dSf_dBeqz_PertQtma, dSt_dBeqz_PertQtma] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertQtma
        %dAbr_dBeqzPertQtma evaluated in x+pert        
        [dAf_dBeqz_PertQtma, dAt_dBeqz_PertQtma] = dAbr_dBeq(dSf_dBeqz_PertQtma, dSt_dBeqz_PertQtma, Sf_PertQtma, St_PertQtma);             
        %2nd Derivatives of Sbr w.r.t. BeqzQtma
        d2Af_dQtmaBeqz(:, k) = (dAf_dBeqz_PertQtma - dAf_dBeqz).' * lam / pert;  %BeqzQtma from, size of [nBeqz, nQtma]
        d2At_dQtmaBeqz(:, k) = (dAt_dBeqz_PertQtma - dAt_dBeqz).' * lam / pert;  %BeqzQtma  to , size of [nBeqz, nQtma]        
    end
    %BeqvQtma num_G46
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
        %dSbr_dBeqvPertQtma evaluated in x+pert
        [dSf_dBeqv_PertQtma, dSt_dBeqv_PertQtma] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertQtma
        %dAbr_dBeqvPertQtma evaluated in x+pert        
        [dAf_dBeqv_PertQtma, dAt_dBeqv_PertQtma] = dAbr_dBeq(dSf_dBeqv_PertQtma, dSt_dBeqv_PertQtma, Sf_PertQtma, St_PertQtma);             
        %2nd Derivatives of Sbr w.r.t. BeqvQtma
        d2Af_dQtmaBeqv(:, k) = (dAf_dBeqv_PertQtma - dAf_dBeqv).' * lam / pert;  %BeqvQtma from, size of [nBeqv, nQtma]
        d2At_dQtmaBeqv(:, k) = (dAt_dBeqv_PertQtma - dAt_dBeqv).' * lam / pert;  %BeqvQtma  to , size of [nBeqv, nQtma]        
    end
    %PfshQtma num_G56
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertQtma, St_PertQtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfshPertQtma evaluated in x+pert
        [dSf_dPfsh_PertQtma, dSt_dPfsh_PertQtma] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertQtma
        %dAbr_dPfshPertQtma evaluated in x+pert        
        [dAf_dPfsh_PertQtma, dAt_dPfsh_PertQtma] = dAbr_dsh(dSf_dPfsh_PertQtma, dSt_dPfsh_PertQtma, Sf_PertQtma, St_PertQtma);             
        %2nd Derivatives of Sbr w.r.t. PfshQtma
        d2Af_dQtmaPfsh(:, k) = (dAf_dPfsh_PertQtma - dAf_dPfsh).' * lam / pert;  %PfshQtma from, size of [nPfsh, nQtma]
        d2At_dQtmaPfsh(:, k) = (dAt_dPfsh_PertQtma - dAt_dPfsh).' * lam / pert;  %PfshQtma  to , size of [nPfsh, nQtma]        
    end
    %Qtma2 num_G66
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertQtma, St_PertQtma] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);       
        %dSbr_dQtmaPertQtma evaluated in x+pert
        [dSf_dQtma_PertQtma, dSt_dQtma_PertQtma] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertQtma
        %dAbr_dQtmaPertQtma evaluated in x+pert        
        [dAf_dQtma_PertQtma, dAt_dQtma_PertQtma] = dAbr_dma(dSf_dQtma_PertQtma, dSt_dQtma_PertQtma, Sf_PertQtma, St_PertQtma);             
        %2nd Derivatives of Abr w.r.t. Qtma2   
        d2Af_dQtma2(:, k) = (dAf_dQtma_PertQtma - dAf_dQtma).' * lam / pert;  %Qtma2 from, size of [nQtma , nQtma] 
        d2At_dQtma2(:, k) = (dAt_dQtma_PertQtma - dAt_dQtma).' * lam / pert;  %Qtma2  to , size of [nQtma , nQtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf16 = sparse(d2Af_dQtmaVa);
num_Hf26 = sparse(d2Af_dQtmaVm);
num_Hf36 = sparse(d2Af_dQtmaBeqz);
num_Hf46 = sparse(d2Af_dQtmaBeqv);
num_Hf56 = sparse(d2Af_dQtmaPfsh);

num_Hf61 = sparse(d2Af_dVaQtma);
num_Hf62 = sparse(d2Af_dVmQtma);
num_Hf63 = sparse(d2Af_dBeqzQtma);
num_Hf64 = sparse(d2Af_dBeqvQtma);
num_Hf65 = sparse(d2Af_dPfshQtma);

num_Hf66 = sparse(d2Af_dQtma2);

num_Ht16 = sparse(d2At_dQtmaVa);
num_Ht26 = sparse(d2At_dQtmaVm);
num_Ht36 = sparse(d2At_dQtmaBeqz);
num_Ht46 = sparse(d2At_dQtmaBeqv);
num_Ht56 = sparse(d2At_dQtmaPfsh);

num_Ht61 = sparse(d2At_dVaQtma);
num_Ht62 = sparse(d2At_dVmQtma);
num_Ht63 = sparse(d2At_dBeqzQtma);
num_Ht64 = sparse(d2At_dBeqvQtma);
num_Ht65 = sparse(d2At_dPfshQtma);

num_Ht66 = sparse(d2At_dQtma2);
