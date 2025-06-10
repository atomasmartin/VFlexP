function [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33,...
          num_Ht13, num_Ht23, num_Ht31, num_Ht32, num_Ht33] = d2Abr_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXBEQZ2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. BeqzVa, BeqzVm, BeqzBeqz and contrariwise (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%   (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [NUM_HF13, NUM_HF23, NUM_HF31, NUM_HF32, NUM_HF33,...
%    NUM_HT13, NUM_HT23, NUM_HT31, NUM_HT32, NUM_HT33] = D2ABR_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 10 matrices containing the 2nd partial derivatives of Abr
%   where:
%   Abr = abs(Sbr).^2
%       = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   Second Derivative w.r.t Beqz 
%       d2Af/dxBeq = (dAfPert - dAf)/pert 
%       d2At/dxBeq = (dAtPert - dAt)/pert
%   thus,
%
%   H13 = HBzVa = d2Abr_dBeqzVa   = d/dBeqz ( (dAbr/dVa  )'*mu )
%   H23 = HBzVm = d2Abr_dBeqzVm   = d/dBeqz ( (dAbr/dVm  )'*mu )
%   H31 = HVaBz = d2Abr_dVaBeqz   = d/dVa   ( (dAbr/dBeqz)'*mu )
%   H32 = HVmBz = d2Abr_dVmBeqz   = d/dVm   ( (dAbr/dBeqz)'*mu )
%   H33 = HBzBz = d2Abr_dBeqz2    = d/dBeqz ( (dAbr/dBeqz)'*mu )
%   
%   H13 = HBzVa = d2Abr_dBeqzVa   =  (  2*real(d2Sbr_dBeqzVa )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dBeqz ))'*mu 
%   H23 = HBzVm = d2Abr_dBeqzVm   =  (  2*real(d2Sbr_dBeqzVm )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dBeqz ))'*mu
%   H31 = HVaBz = d2Abr_dVaBeqz   =  (  2*real(d2Sbr_dVaBeqz )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dVa   ))'*mu
%   H32 = HVmBz = d2Abr_dVmBeqz   =  (  2*real(d2Sbr_dVmBeqz )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dVm   ))'*mu
%   H33 = HBzBz = d2Abr_dBeqz2    =  (  2*real(d2Sbr_dBeqz2  )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dBeqz ))'*mu
%
%   Example:
%   [NUM_HF13, NUM_HF23, NUM_HF31, NUM_HF32, NUM_HF33,...
%    NUM_HT13, NUM_HT23, NUM_HT31, NUM_HT32, NUM_HT33] = D2ABR_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2.
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

if vcart
    error('d2Abr_dxBeqz2Pert: Derivatives of Flow Limit equations w.r.t Beq using Finite Differences Method in cartasian have not been coded yet')    

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
    
    %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    
    %Dimensionalize (Allocate for computational speed) 
    d2Af_dBeqzVa = sparse( zeros(nb,   nBeqz) );
    d2Af_dBeqzVm = sparse( zeros(nb,   nBeqz) );

    d2Af_dVaBeqz = sparse( zeros(nBeqz,   nb) );
    d2Af_dVmBeqz = sparse( zeros(nBeqz,   nb) ); 

    d2Af_dBeqz2  = sparse( zeros(nBeqz,nBeqz) );
    
    d2At_dBeqzVa = sparse( zeros(nb,   nBeqz) );
    d2At_dBeqzVm = sparse( zeros(nb,   nBeqz) );

    d2At_dVaBeqz = sparse( zeros(nBeqz,   nb) );
    d2At_dVmBeqz = sparse( zeros(nBeqz,   nb) );     
    
    d2At_dBeqz2  = sparse( zeros(nBeqz,nBeqz) );  
    
    %BeqzVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dBeqz_PertVa, dSt_dBeqz_PertVa] = dSbr_dBeq(branch, V1p, 1, vcart); %Sbr_PertVa
        [dAf_dBeqz_PertVa, dAt_dBeqz_PertVa] = dAbr_dBeq(dSf_dBeqz_PertVa, dSt_dBeqz_PertVa, Sf_PertVa, St_PertVa); %dAbr_dBeqzPertVa
        d2Af_dVaBeqz(:, k) = (dAf_dBeqz_PertVa - dAf_dBeqz).' * lam / pert; %BeqzVa From
        d2At_dVaBeqz(:, k) = (dAt_dBeqz_PertVa - dAt_dBeqz).' * lam / pert; %BeqzVa To
    end
    %BeqzVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dBeqz_PertVm, dSt_dBeqz_PertVm] = dSbr_dBeq(branch, V2p, 1, vcart); %Sbr_PertVm
        [dAf_dBeqz_PertVm, dAt_dBeqz_PertVm] = dAbr_dBeq(dSf_dBeqz_PertVm, dSt_dBeqz_PertVm, Sf_PertVm, St_PertVm); %dAbr_dBeqzPertVm
        d2Af_dVmBeqz(:, k) = (dAf_dBeqz_PertVm - dAf_dBeqz).' * lam / pert; %BeqzVm From
        d2At_dVmBeqz(:, k) = (dAt_dBeqz_PertVm - dAt_dBeqz).' * lam / pert; %BeqzVm To
    end
    
    %VxBeqz
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVaPertBeqz evaluated in x+pert        
        [dSf_dV1_PertBeqz, dSf_dV2_PertBeqz, dSt_dV1_PertBeqz, dSt_dV2_PertBeqz, Sf_PertBeqz, St_PertBeqz] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVaPertBeqz evaluated in x+pert
        [dAf_dV1_PertBeqz, dAf_dV2_PertBeqz, dAt_dV1_PertBeqz, dAt_dV2_PertBeqz] = dAbr_dV(dSf_dV1_PertBeqz, dSf_dV2_PertBeqz, dSt_dV1_PertBeqz, dSt_dV2_PertBeqz, Sf_PertBeqz, St_PertBeqz);
        %2nd Derivatives of Abr w.r.t. BeqzVx
        d2Af_dBeqzVa(:, k) = (dAf_dV1_PertBeqz - dAf_dV1).' * lam / pert;  %VaBeqz from, size of [nb, nBeqz] 
        d2Af_dBeqzVm(:, k) = (dAf_dV2_PertBeqz - dAf_dV2).' * lam / pert;  %VmBeqz from, size of [nb, nBeqz]
        d2At_dBeqzVa(:, k) = (dAt_dV1_PertBeqz - dAt_dV1).' * lam / pert;  %VaBeqz  to , size of [nb, nBeqz] 
        d2At_dBeqzVm(:, k) = (dAt_dV2_PertBeqz - dAt_dV2).' * lam / pert;  %VmBeqz  to , size of [nb, nBeqz]        
    end
    %BeqzBeqz
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
        %dSbr_dBeqzPertBeqz evaluated in x+pert
        [dSf_dBeqz_PertBeqz, dSt_dBeqz_PertBeqz] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertBeqz
        %dAbr_dBeqzPertBeqz evaluated in x+pert        
        [dAf_dBeqz_PertBeqz, dAt_dBeqz_PertBeqz] = dAbr_dBeq(dSf_dBeqz_PertBeqz, dSt_dBeqz_PertBeqz, Sf_PertBeqz, St_PertBeqz);       
        %2nd Derivatives of Abr w.r.t. Beqz2   
        d2Af_dBeqz2(:, k) = (dAf_dBeqz_PertBeqz - dAf_dBeqz).' * lam / pert;  %BeqzBeqz from, size of [nBeqz , nBeqz] 
        d2At_dBeqz2(:, k) = (dAt_dBeqz_PertBeqz - dAt_dBeqz).' * lam / pert;  %BeqzBeqz  to , size of [nBeqz , nBeqz]
    end
end

num_Hf13 = sparse(d2Af_dBeqzVa);
num_Hf23 = sparse(d2Af_dBeqzVm);

num_Hf31 = sparse(d2Af_dVaBeqz);
num_Hf32 = sparse(d2Af_dVmBeqz);

num_Hf33 = sparse(d2Af_dBeqz2);

num_Ht13 = sparse(d2At_dBeqzVa);
num_Ht23 = sparse(d2At_dBeqzVm);

num_Ht31 = sparse(d2At_dVaBeqz);
num_Ht32 = sparse(d2At_dVmBeqz);

num_Ht33 = sparse(d2At_dBeqz2);
