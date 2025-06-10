function [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44,...
          num_Ht14, num_Ht24, num_Ht34, num_Ht41, num_Ht42, num_Ht43, num_Ht44] = d2Abr_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXBEQV2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. BeqvVa, BeqvVm, BeqvBeqz, BeqvBeqv and contrariwise (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%   (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [NUM_HF14, NUM_HF24, NUM_HF34, NUM_HF41, NUM_HF42, NUM_HF43, NUM_HF44,...
%    NUM_HT14, NUM_HT24, NUM_HT34, NUM_HT41, NUM_HT42, NUM_HT43, NUM_HT44] = D2ABR_DXBEQV2(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 14 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H14 = HBvVa = d2Abr_dBeqvVa   = d/dBeqv  ( (dAbr/dVa  )'*mu )
%   H24 = HBvVm = d2Abr_dBeqvVm   = d/dBeqv  ( (dAbr/dVm  )'*mu )
%   H34 = HBvBz = d2Abr_dBeqvBeqz = d/dBeqv  ( (dAbr/dBeqz)'*mu )
%   H41 = HVaBv = d2Abr_dVaBeqv   = d/dVa    ( (dAbr/dBeqv)'*mu )
%   H42 = HVmBv = d2Abr_dVmBeqv   = d/dVm    ( (dAbr/dBeqv)'*mu )
%   H43 = HBzBv = d2Abr_dBeqzBeqv = d/dBeqz  ( (dAbr/dBeqv)'*mu )
%   H44 = HBvBv = d2Abr_dBeqvBeqv = d/dBeqv  ( (dAbr/dBeqv)'*mu )
%   
%   H14 = HBvVa = d2Abr_dBeqvVa   =  (  2*real(d2Sbr_dBeqvVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dBeqv ))'*mu 
%   H24 = HBvVm = d2Abr_dBeqvVm   =  (  2*real(d2Sbr_dBeqvVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dBeqv ))'*mu
%   H34 = HBvBz = d2Abr_dBeqvBeqz =  (  2*real(d2Sbr_dBeqvBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dBeqv ))'*mu
%   H41 = HVaBv = d2Abr_dVaBeqv   =  (  2*real(d2Sbr_dVaBeqv  )*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dVa   ))'*mu
%   H42 = HVmBv = d2Abr_dVmBeqv   =  (  2*real(d2Sbr_dVmBeqv  )*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dVm   ))'*mu
%   H43 = HBzBv = d2Abr_dBeqzBeqv =  (  2*real(d2Sbr_dBeqzBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dBeqz ))'*mu
%   H44 = HBvBv = d2Abr_dBeqvBeqv =  (  2*real(d2Sbr_dBeqvBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dBeqv ))'*mu
%
%   Example:
%   [NUM_HF14, NUM_HF24, NUM_HF34, NUM_HF41, NUM_HF42, NUM_HF43, NUM_HF44,...
%    NUM_HT14, NUM_HT24, NUM_HT34, NUM_HT41, NUM_HT42, NUM_HT43, NUM_HT44] = D2ABR_DXBEQV2(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)%
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXBEQZ2.
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
%[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Calculation of derivatives
if vcart
    error('d2Abr_dxBeqv2Pert: Derivatives of Flow Limit equations w.r.t Beq using Finite Differences Method in cartasian have not been coded yet')    

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
    
    %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St); 
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
      
    %Selector of active Beq 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    
    %Dimensionalize (Allocate for computational speed) 
    d2Af_dBeqvVa   = sparse( zeros(nb   ,nBeqv) );
    d2Af_dBeqvVm   = sparse( zeros(nb   ,nBeqv) );
    d2Af_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    d2Af_dVaBeqv   = sparse( zeros(nBeqv,   nb) );
    d2Af_dVmBeqv   = sparse( zeros(nBeqv,   nb) ); 
    d2Af_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );
    
    d2Af_dBeqv2    = sparse( zeros(nBeqv,nBeqv) );
    
    d2At_dBeqvVa   = sparse( zeros(nb   ,nBeqv) );
    d2At_dBeqvVm   = sparse( zeros(nb   ,nBeqv) );
    d2At_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    d2At_dVaBeqv   = sparse( zeros(nBeqv,   nb) );
    d2At_dVmBeqv   = sparse( zeros(nBeqv,   nb) ); 
    d2At_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );
    
    d2At_dBeqv2    = sparse( zeros(nBeqv,nBeqv) ); 

    %BeqvVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dBeqv_PertVa, dSt_dBeqv_PertVa] = dSbr_dBeq(branch, V1p, 2, vcart); %dSbr_dBeqvPertVa
        [dAf_dBeqv_PertVa, dAt_dBeqv_PertVa] = dAbr_dBeq(dSf_dBeqv_PertVa, dSt_dBeqv_PertVa, Sf_PertVa, St_PertVa); %dAbr_dBeqvPertVa       
        d2Af_dVaBeqv(:, k) = (dAf_dBeqv_PertVa - dAf_dBeqv).' * lam / pert; %BeqvVa From
        d2At_dVaBeqv(:, k) = (dAt_dBeqv_PertVa - dAt_dBeqv).' * lam / pert; %BeqvVa To
   end
    %BeqvVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dBeqv_PertVm, dSt_dBeqv_PertVm] = dSbr_dBeq(branch, V2p, 2, vcart); %dSbr_dBeqvPertVm
        [dAf_dBeqv_PertVm, dAt_dBeqv_PertVm] = dAbr_dBeq(dSf_dBeqv_PertVm, dSt_dBeqv_PertVm, Sf_PertVm, St_PertVm); %dAbr_dBeqvPertVm       
        d2Af_dVmBeqv(:, k) = (dAf_dBeqv_PertVm - dAf_dBeqv).' * lam / pert; %BeqvVm From
        d2At_dVmBeqv(:, k) = (dAt_dBeqv_PertVm - dAt_dBeqv).' * lam / pert; %BeqvVm To
   end
    %BeqvBeqz
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
        %dSbr_dBeqvPertBeqz evaluated in x+pert
        [dSf_dBeqv_PertBeqz, dSt_dBeqv_PertBeqz] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertBeqz
        %dAbr_dBeqvPertBeqz evaluated in x+pert        
        [dAf_dBeqv_PertBeqz, dAt_dBeqv_PertBeqz] = dAbr_dBeq(dSf_dBeqv_PertBeqz, dSt_dBeqv_PertBeqz, Sf_PertBeqz, St_PertBeqz); %dAbr_dBeqvPertBeqz              
        %2nd Derivatives of Abr w.r.t. BeqvBeqz
        d2Af_dBeqzBeqv(:, k) = (dAf_dBeqv_PertBeqz - dAf_dBeqv).' * lam / pert; %BeqvBeqz From
        d2At_dBeqzBeqv(:, k) = (dAt_dBeqv_PertBeqz - dAt_dBeqv).' * lam / pert; %BeqvBeqz To 
    end
    %VxBeqv
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertBeqv evaluated in x+pert        
        [dSf_dV1_PertBeqv, dSf_dV2_PertBeqv, dSt_dV1_PertBeqv, dSt_dV2_PertBeqv, Sf_PertBeqv, St_PertBeqv] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVxPertBeqv evaluated in x+pert
        [dAf_dV1_PertBeqv, dAf_dV2_PertBeqv, dAt_dV1_PertBeqv, dAt_dV2_PertBeqv] = dAbr_dV(dSf_dV1_PertBeqv, dSf_dV2_PertBeqv, dSt_dV1_PertBeqv, dSt_dV2_PertBeqv, Sf_PertBeqv, St_PertBeqv);
        %2nd Derivatives of Abr w.r.t. BeqvVx
        d2Af_dBeqvVa(:, k) = (dAf_dV1_PertBeqv - dAf_dV1).' * lam / pert;  %VaBeqv from, size of [nb, nBeqv] 
        d2Af_dBeqvVm(:, k) = (dAf_dV2_PertBeqv - dAf_dV2).' * lam / pert;  %VmBeqv from, size of [nb, nBeqv]
        d2At_dBeqvVa(:, k) = (dAt_dV1_PertBeqv - dAt_dV1).' * lam / pert;  %VaBeqv  to , size of [nb, nBeqv] 
        d2At_dBeqvVm(:, k) = (dAt_dV2_PertBeqv - dAt_dV2).' * lam / pert;  %VmBeqv  to , size of [nb, nBeqv]        
    end
    %BeqzBeqv
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
        %dSbr_dBeqzPertBeqv evaluated in x+pert
        [dSf_dBeqz_PertBeqv, dSt_dBeqz_PertBeqv] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertBeqv
        %dAbr_dBeqzPertBeqv evaluated in x+pert        
        [dAf_dBeqz_PertBeqv, dAt_dBeqz_PertBeqv] = dAbr_dBeq(dSf_dBeqz_PertBeqv, dSt_dBeqz_PertBeqv, Sf_PertBeqv, St_PertBeqv);             
        %2nd Derivatives of Sbr w.r.t. BeqzBeqv
        d2Af_dBeqvBeqz(:, k) = (dAf_dBeqz_PertBeqv - dAf_dBeqz).' * lam / pert;  %BeqzBeqv from, size of [nBeqz, nBeqv]
        d2At_dBeqvBeqz(:, k) = (dAt_dBeqz_PertBeqv - dAt_dBeqz).' * lam / pert;  %BeqzBeqv  to , size of [nBeqz, nBeqv]        
    end
    %BeqvBeqv
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
        %dSbr_dBeqvPertBeqv evaluated in x+pert
        [dSf_dBeqv_PertBeqv, dSt_dBeqv_PertBeqv] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertBeqv
        %dAbr_dBeqvPertBeqv evaluated in x+pert        
        [dAf_dBeqv_PertBeqv, dAt_dBeqv_PertBeqv] = dAbr_dBeq(dSf_dBeqv_PertBeqv, dSt_dBeqv_PertBeqv, Sf_PertBeqv, St_PertBeqv);       
        %2nd Derivatives of Abr w.r.t. Beqv2   
        d2Af_dBeqv2(:, k) = (dAf_dBeqv_PertBeqv - dAf_dBeqv).' * lam / pert;  %BeqvBeqv from, size of [nBeqv , nBeqv] 
        d2At_dBeqv2(:, k) = (dAt_dBeqv_PertBeqv - dAt_dBeqv).' * lam / pert;  %BeqvBeqv  to , size of [nBeqv , nBeqv]
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf14 = sparse(d2Af_dBeqvVa);
num_Hf24 = sparse(d2Af_dBeqvVm);
num_Hf34 = sparse(d2Af_dBeqvBeqz);

num_Hf41 = sparse(d2Af_dVaBeqv);
num_Hf42 = sparse(d2Af_dVmBeqv);
num_Hf43 = sparse(d2Af_dBeqzBeqv);

num_Hf44 = sparse(d2Af_dBeqv2);

num_Ht14 = sparse(d2At_dBeqvVa);
num_Ht24 = sparse(d2At_dBeqvVm);
num_Ht34 = sparse(d2At_dBeqvBeqz);

num_Ht41 = sparse(d2At_dVaBeqv);
num_Ht42 = sparse(d2At_dVmBeqv);
num_Ht43 = sparse(d2At_dBeqzBeqv);
num_Ht44 = sparse(d2At_dBeqv2);

