function [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55,...
          num_Ht15, num_Ht25, num_Ht35, num_Ht45, num_Ht51, num_Ht52, num_Ht53, num_Ht54, num_Ht55] = d2Abr_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2ABR_DXSH2PERT  Computes 2nd derivatives of |branch flow|^2 w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%  [NUM_HF15, NUM_HF25, NUM_HF35, NUM_HF45, NUM_HF51, NUM_HF52, NUM_HF53, NUM_HF54, NUM_HF55,...
%   NUM_HT15, NUM_HT25, NUM_HT35, NUM_HT45, NUM_HT51, NUM_HT52, NUM_HT53, NUM_HT54, NUM_HT55] = D2ABR_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 18 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H15 = HShVa = d2Abr_dShVa     = d/dsh   ( (dAbr/dVa  )'*mu )
%   H25 = HShVm = d2Abr_dShVm     = d/dsh   ( (dAbr/dVm  )'*mu )
%   H35 = HShBz = d2Abr_dShBeqz   = d/dsh   ( (dAbr/dBeqz)'*mu )
%   H45 = HShBv = d2Abr_dShBeqv   = d/dsh   ( (dAbr/dBeqv)'*mu )
%   H51 = HVaSh = d2Abr_dVaSh     = d/dVa   ( (dAbr/dsh  )'*mu )
%   H52 = HVmSh = d2Abr_dVmSh     = d/dVm   ( (dAbr/dsh  )'*mu )
%   H53 = HBzSh = d2Abr_dBeqzSh   = d/dBeqz ( (dAbr/dsh  )'*mu )
%   H54 = HBvSh = d2Abr_dBeqvSh   = d/dBeqv ( (dAbr/dsh  )'*mu )
%   H55 = HShSh = d2Abr_dShSh     = d/dsh   ( (dAbr/dsh  )'*mu )
%   
%   H15 = HShVa = d2Abr_dShVa     =  (  2*real(d2Sbr_dShVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dSh   ))'*mu 
%   H25 = HShVm = d2Abr_dShVm     =  (  2*real(d2Sbr_dShVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dSh   ))'*mu
%   H35 = HShBz = d2Abr_dShBeqz   =  (  2*real(d2Sbr_dShBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dSh   ))'*mu
%   H45 = HShBv = d2Abr_dShBeqv   =  (  2*real(d2Sbr_dShBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dSh   ))'*mu
%   H51 = HVaSh = d2Abr_dVaSh     =  (  2*real(d2Sbr_dVaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dVa   ))'*mu
%   H52 = HVmSh = d2Abr_dVmSh     =  (  2*real(d2Sbr_dVmSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dVm   ))'*mu
%   H53 = HBzSh = d2Abr_dBeqzSh   =  (  2*real(d2Sbr_dBeqzSh)*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dBeqz ))'*mu
%   H54 = HBvSh = d2Abr_dBeqvSh   =  (  2*real(d2Sbr_dBeqvSh)*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dBeqv ))'*mu
%   H55 = HShSh = d2Abr_dShSh     =  (  2*real(d2Sbr_dShSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dSh   ))'*mu 
%
%   Example:
%  [NUM_HF15, NUM_HF25, NUM_HF35, NUM_HF45, NUM_HF51, NUM_HF52, NUM_HF53, NUM_HF54, NUM_HF55,...
%   NUM_HT15, NUM_HT25, NUM_HT35, NUM_HT45, NUM_HT51, NUM_HT52, NUM_HT53, NUM_HT54, NUM_HT55] = D2ABR_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXBEQV2.
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
%[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift

%% Calculation of derivatives
if vcart
    error('d2Abr_dxsh2Pert: Derivatives of Flow Limit equations w.r.t Theta_sh using Finite Differences Method in cartasian have not been coded yet')    

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
    
     %Abr 1st Derivatives    
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [dAf_dBeqz, dAt_dBeqz] = dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St); 
    [dAf_dPfsh, dAt_dPfsh] = dAbr_dsh(dSf_dPfsh, dSt_dPfsh, Sf, St); 
        
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
   
    %Dimensionalize (Allocate for computational speed)  
    d2Af_dPfshVa   = sparse( zeros(nb,   nPfsh) ); 
    d2Af_dPfshVm   = sparse( zeros(nb,   nPfsh) );
    d2Af_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2Af_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );

    d2Af_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2Af_dVmPfsh   = sparse( zeros(nPfsh,nb   ) );
    d2Af_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2Af_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
    
    d2Af_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );

    d2At_dPfshVa   = sparse( zeros(nb,   nPfsh) ); 
    d2At_dPfshVm   = sparse( zeros(nb,   nPfsh) );
    d2At_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2At_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );

    d2At_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2At_dVmPfsh   = sparse( zeros(nPfsh,nb   ) );
    d2At_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2At_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
   
    d2At_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );    
    
    %PfShVa num_G51
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [Sf_PertVa, St_PertVa] = SbrFlows(branch, Yf, Yt, V1p);
        [dSf_dPfsh_PertVa, dSt_dPfsh_PertVa] = dSbr_dsh(branch, V1p, 1, vcart); %dSbr_dPfshPertVa
        [dAf_dPfsh_PertVa, dAt_dPfsh_PertVa] = dAbr_dsh(dSf_dPfsh_PertVa, dSt_dPfsh_PertVa, Sf_PertVa, St_PertVa); %dAbr_dPfshPertVa       
        d2Af_dVaPfsh(:, k) = (dAf_dPfsh_PertVa - dAf_dPfsh).' * lam / pert; %shVa from, size of [nPfsh, nb]
        d2At_dVaPfsh(:, k) = (dAt_dPfsh_PertVa - dAt_dPfsh).' * lam / pert; %shVa  to , size of [nPfsh, nb]
    end
    %PfshVm num_G52
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [Sf_PertVm, St_PertVm] = SbrFlows(branch, Yf, Yt, V2p);
        [dSf_dPfsh_PertVm, dSt_dPfsh_PertVm] = dSbr_dsh(branch, V2p, 1, vcart); %dSbr_dPfshPertVm
        [dAf_dPfsh_PertVm, dAt_dPfsh_PertVm] = dAbr_dsh(dSf_dPfsh_PertVm, dSt_dPfsh_PertVm, Sf_PertVm, St_PertVm); %dAbr_dPfshPertVm       
        d2Af_dVmPfsh(:, k) = (dAf_dPfsh_PertVm - dAf_dPfsh).' * lam / pert; %shVm from, size of [nPfsh, nb]
        d2At_dVmPfsh(:, k) = (dAt_dPfsh_PertVm - dAt_dPfsh).' * lam / pert; %shVm  to , size of [nPfsh, nb]
    end
    %PfshBeqz num_G53
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
        %dSbr_dPfshPertBeqz evaluated in x+pert
        [dSf_dPfsh_PertBeqz, dSt_dPfsh_PertBeqz] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertBeqz
        %dAbr_dPfshPertBeqz evaluated in x+pert        
        [dAf_dPfsh_PertBeqz, dAt_dPfsh_PertBeqz] = dAbr_dsh(dSf_dPfsh_PertBeqz, dSt_dPfsh_PertBeqz, Sf_PertBeqz, St_PertBeqz); %dAbr_dPfshPertBeqz              
        %2nd Derivatives of Abr w.r.t. PfshBeqz
        d2Af_dBeqzPfsh(:, k) = (dAf_dPfsh_PertBeqz - dAf_dPfsh).' * lam / pert; %PfshBeqz From
        d2At_dBeqzPfsh(:, k) = (dAt_dPfsh_PertBeqz - dAt_dPfsh).' * lam / pert; %PfshBeqz To 
    end
    %PfshBeqv num_G54
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
        %dSbr_dPfshPertBeqv evaluated in x+pert
        [dSf_dPfsh_PertBeqv, dSt_dPfsh_PertBeqv] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertBeqv
        %dAbr_dPfshPertBeqv evaluated in x+pert        
        [dAf_dPfsh_PertBeqv, dAt_dPfsh_PertBeqv] = dAbr_dsh(dSf_dPfsh_PertBeqv, dSt_dPfsh_PertBeqv, Sf_PertBeqv, St_PertBeqv); %dAbr_dPfshPertBeqv              
        %2nd Derivatives of Abr w.r.t. PfshBeqv
        d2Af_dBeqvPfsh(:, k) = (dAf_dPfsh_PertBeqv - dAf_dPfsh).' * lam / pert; %PfshBeqv From
        d2At_dBeqvPfsh(:, k) = (dAt_dPfsh_PertBeqv - dAt_dPfsh).' * lam / pert; %PfshBeqv To 
    end
    %VxPfsh num_G15 num_G25
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertPfsh evaluated in x+pert        
        [dSf_dV1_PertPfsh, dSf_dV2_PertPfsh, dSt_dV1_PertPfsh, dSt_dV2_PertPfsh, Sf_PertPfsh, St_PertPfsh] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %dAbr_dVxPertPfsh evaluated in x+pert
        [dAf_dV1_PertPfsh, dAf_dV2_PertPfsh, dAt_dV1_PertPfsh, dAt_dV2_PertPfsh] = dAbr_dV(dSf_dV1_PertPfsh, dSf_dV2_PertPfsh, dSt_dV1_PertPfsh, dSt_dV2_PertPfsh, Sf_PertPfsh, St_PertPfsh);
        %2nd Derivatives of Abr w.r.t. PfshVx
        d2Af_dPfshVa(:, k) = (dAf_dV1_PertPfsh - dAf_dV1).' * lam / pert;  %VaPfsh from, size of [nb, nPfsh] 
        d2Af_dPfshVm(:, k) = (dAf_dV2_PertPfsh - dAf_dV2).' * lam / pert;  %VmPfsh from, size of [nb, nPfsh]
        d2At_dPfshVa(:, k) = (dAt_dV1_PertPfsh - dAt_dV1).' * lam / pert;  %VaPfsh  to , size of [nb, nPfsh] 
        d2At_dPfshVm(:, k) = (dAt_dV2_PertPfsh - dAt_dV2).' * lam / pert;  %VmPfsh  to , size of [nb, nPfsh]        
    end
    %BeqzPfsh num_G35
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
        %dSbr_dBeqzPertPfsh evaluated in x+pert
        [dSf_dBeqz_PertPfsh, dSt_dBeqz_PertPfsh] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertPfsh
        %dAbr_dBeqzPertPfsh evaluated in x+pert        
        [dAf_dBeqz_PertPfsh, dAt_dBeqz_PertPfsh] = dAbr_dBeq(dSf_dBeqz_PertPfsh, dSt_dBeqz_PertPfsh, Sf_PertPfsh, St_PertPfsh);             
        %2nd Derivatives of Sbr w.r.t. BeqzPfsh
        d2Af_dPfshBeqz(:, k) = (dAf_dBeqz_PertPfsh - dAf_dBeqz).' * lam / pert;  %BeqzPfsh from, size of [nBeqz, nPfsh]
        d2At_dPfshBeqz(:, k) = (dAt_dBeqz_PertPfsh - dAt_dBeqz).' * lam / pert;  %BeqzPfsh  to , size of [nBeqz, nPfsh]        
    end
    %BeqvPfsh num_G45
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
        %dSbr_dBeqvPertPfsh evaluated in x+pert
        [dSf_dBeqv_PertPfsh, dSt_dBeqv_PertPfsh] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertPfsh
        %dAbr_dBeqvPertPfsh evaluated in x+pert        
        [dAf_dBeqv_PertPfsh, dAt_dBeqv_PertPfsh] = dAbr_dBeq(dSf_dBeqv_PertPfsh, dSt_dBeqv_PertPfsh, Sf_PertPfsh, St_PertPfsh);   
        %2nd Derivatives of Sbr w.r.t. BeqvPfsh
        d2Af_dPfshBeqv(:, k) = (dAf_dBeqv_PertPfsh - dAf_dBeqv).' * lam / pert;  %BeqvPfsh from, size of [nBeqv, nPfsh]
        d2At_dPfshBeqv(:, k) = (dAt_dBeqv_PertPfsh - dAt_dBeqv).' * lam / pert;  %BeqvPfsh  to , size of [nBeqv, nPfsh]        
    end
    %PfshPfsh num_G55
    for k=1:nPfsh
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShSel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %Sbr evaluated in x+pert
        [Sf_PertPfsh, St_PertPfsh] = SbrFlows(branch_Pert, Yf_Pert, Yt_Pert, V);
        %dSbr_dPfshPertPfsh evaluated in x+pert
        [dSf_dPfsh_PertPfsh, dSt_dPfsh_PertPfsh] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertPfsh
        %dAbr_dPfshPertPfsh evaluated in x+pert        
        [dAf_dPfsh_PertPfsh, dAt_dPfsh_PertPfsh] = dAbr_dsh(dSf_dPfsh_PertPfsh, dSt_dPfsh_PertPfsh, Sf_PertPfsh, St_PertPfsh);       
        %2nd Derivatives of Abr w.r.t. Pfsh2   
        d2Af_dPfsh2(:, k) = (dAf_dPfsh_PertPfsh - dAf_dPfsh).' * lam / pert;  %PfshPfsh from, size of [nPfsh , nPfsh] 
        d2At_dPfsh2(:, k) = (dAt_dPfsh_PertPfsh - dAt_dPfsh).' * lam / pert;  %PfshPfsh  to , size of [nPfsh , nPfsh]
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf15 = sparse(d2Af_dPfshVa);
num_Hf25 = sparse(d2Af_dPfshVm);
num_Hf35 = sparse(d2Af_dPfshBeqz);
num_Hf45 = sparse(d2Af_dPfshBeqv);

num_Hf51 = sparse(d2Af_dVaPfsh);
num_Hf52 = sparse(d2Af_dVmPfsh);
num_Hf53 = sparse(d2Af_dBeqzPfsh);
num_Hf54 = sparse(d2Af_dBeqvPfsh);

num_Hf55 = sparse(d2Af_dPfsh2);

num_Ht15 = sparse(d2At_dPfshVa);
num_Ht25 = sparse(d2At_dPfshVm);
num_Ht35 = sparse(d2At_dPfshBeqz);
num_Ht45 = sparse(d2At_dPfshBeqv);

num_Ht51 = sparse(d2At_dVaPfsh);
num_Ht52 = sparse(d2At_dVmPfsh);
num_Ht53 = sparse(d2At_dBeqzPfsh);
num_Ht54 = sparse(d2At_dBeqvPfsh);

num_Ht55 = sparse(d2At_dPfsh2);
