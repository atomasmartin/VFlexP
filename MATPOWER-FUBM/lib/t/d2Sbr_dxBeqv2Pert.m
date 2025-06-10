function [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44,...
          num_Ht14, num_Ht24, num_Ht34, num_Ht41, num_Ht42, num_Ht43, num_Ht44] = d2Sbr_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBR_DXBEQV2PERT  Computes 2nd derivatives of complex brch power flow w.r.t. BeqvVa, BeqvVm, BeqvBeqz, VaBeqv, VmBeqv, BeqzBeqv, BeqvBeqv (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv,...
%    HtBvVa, HtBvVm, HtBvBz, HtVaBv, HtVmBv, HtBzBv, HtBvBv] = D2SBR_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 14 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv,...
%    HtBvVa, HtBvVm, HtBvBz, HtVaBv, HtVmBv, HtBzBv, HtBvBv] = D2SBR_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv,...
%    HtBvVa, HtBvVm, HtBvBz, HtVaBv, HtVmBv, HtBzBv, HtBvBv] = D2SBR_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%       Here the output matrices correspond to:
%           HfBvVa = d/dBeqv (dSf_dVa.'   * mu)
%           HfBvVm = d/dBeqv (dSf_dVm.'   * mu)
%           HfBvBz = d/dBeqv (dSf_dBeqz.' * mu)
%           HfVaBv = d/dVa   (dSf_dBeqv.' * mu)
%           HfVmBv = d/dVm   (dSf_dBeqv.' * mu)
%           HfBzBv = d/dBeqz (dSf_dBeqv.' * mu)
%           HfBvBv = d/dBeqv (dSf_dBeqv.' * mu)
%
%           HtBvVa = d/dBeqv (dSt_dVa.'   * mu)
%           HtBvVm = d/dBeqv (dSt_dVm.'   * mu)
%           HtBvBz = d/dBeqv (dSt_dBeqz.' * mu)
%           HtVaBv = d/dVa   (dSt_dBeqv.' * mu)
%           HtVmBv = d/dVm   (dSt_dBeqv.' * mu)
%           HtBzBv = d/dBeqz (dSt_dBeqv.' * mu)
%           HtBvBv = d/dBeqv (dSt_dBeqv.' * mu)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf
%   [TNX]  A. Alvarez-Bustos, "AC/DC Power Flows and their Derivatives using
%          FUBM Complex Matrix Notation" MATPOWER Technical Note x, Month 20XX.
%             http://www.pserc.cornell.edu/matpower/
%                                           TNX-OPF-Derivatives-FUBM.pdf   %%AAB- Technical note to be written
                                          
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
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
    error('d2Sbr_dxBeqv2Pert: Derivatives of Flow Limit equations w.r.t Beq using Finite Differences Method in cartasian have not been coded yet')    

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
     
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
      
    %Selector of active Beq 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqvAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed) 
    d2Sf_dBeqvVa   = sparse( zeros(nb   ,nBeqv) );
    d2Sf_dBeqvVm   = sparse( zeros(nb   ,nBeqv) );
    d2Sf_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    
    d2Sf_dVaBeqv   = sparse( zeros(nBeqv,   nb) );
    d2Sf_dVmBeqv   = sparse( zeros(nBeqv,   nb) ); 
    d2Sf_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );

    d2Sf_dBeqv2    = sparse( zeros(nBeqv,nBeqv) );
    
    d2St_dBeqvVa   = sparse( zeros(nb   ,nBeqv) );
    d2St_dBeqvVm   = sparse( zeros(nb   ,nBeqv) );
    d2St_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    
    d2St_dVaBeqv   = sparse( zeros(nBeqv,   nb) );
    d2St_dVmBeqv   = sparse( zeros(nBeqv,   nb) ); 
    d2St_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );

    d2St_dBeqv2    = sparse( zeros(nBeqv,nBeqv) ); 

    %BeqvVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dBeqv_PertVa, dSt_dBeqv_PertVa] = dSbr_dBeq(branch, V1p, 2, vcart); %dSbr_dBeqvPertVa
        d2Sf_dVaBeqv(:, k) = (dSf_dBeqv_PertVa - dSf_dBeqv).' * lam / pert; %BeqvVa From
        d2St_dVaBeqv(:, k) = (dSt_dBeqv_PertVa - dSt_dBeqv).' * lam / pert; %BeqvVa To
    end
    %BeqvVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dBeqv_PertVm, dSt_dBeqv_PertVm] = dSbr_dBeq(branch, V2p, 2, vcart); %dSbr_dBeqvPertVm
        d2Sf_dVmBeqv(:, k) = (dSf_dBeqv_PertVm - dSf_dBeqv).' * lam / pert; %BeqvVm From
        d2St_dVmBeqv(:, k) = (dSt_dBeqv_PertVm - dSt_dBeqv).' * lam / pert; %BeqvVm To  
    end
    %BeqvBeqz
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqvPertBeqz evaluated in x+pert
        [dSf_dBeqv_PertBeqz, dSt_dBeqv_PertBeqz] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertBeqz
        %2nd Derivatives of Sbr w.r.t. BeqvBeqz
        d2Sf_dBeqzBeqv(:, k) = (dSf_dBeqv_PertBeqz - dSf_dBeqv).' * lam / pert; %BeqvBeqz From
        d2St_dBeqzBeqv(:, k) = (dSt_dBeqv_PertBeqz - dSt_dBeqv).' * lam / pert; %BeqvBeqz To 
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
        %dSbus_dVaPertBeqv evaluated in x+pert
        [dSf_dV1_PertBeqv, dSf_dV2_PertBeqv, dSt_dV1_PertBeqv, dSt_dV2_PertBeqv, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. BeqvVx
        d2Sf_dBeqvVa(:, k) = (dSf_dV1_PertBeqv - dSf_dV1).' * lam / pert;  %VaBeqv from, size of [nb, nBeqv] 
        d2Sf_dBeqvVm(:, k) = (dSf_dV2_PertBeqv - dSf_dV2).' * lam / pert;  %VmBeqv from, size of [nb, nBeqv] 
        d2St_dBeqvVa(:, k) = (dSt_dV1_PertBeqv - dSt_dV1).' * lam / pert;  %VaBeqv  to , size of [nb, nBeqv] 
        d2St_dBeqvVm(:, k) = (dSt_dV2_PertBeqv - dSt_dV2).' * lam / pert;  %VmBeqv  to , size of [nb, nBeqv]         
    end
    %BeqzBeqv
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertBeqv evaluated in x+pert
        [dSf_dBeqz_PertBeqv, dSt_dBeqz_PertBeqv] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertBeqv
        %2nd Derivatives of Sbr w.r.t. BeqzBeqv
        d2Sf_dBeqvBeqz(:, k) = (dSf_dBeqz_PertBeqv - dSf_dBeqz).' * lam / pert;  %BeqzBeqv from, size of [nBeqz, nBeqv]
        d2St_dBeqvBeqz(:, k) = (dSt_dBeqz_PertBeqv - dSt_dBeqz).' * lam / pert;  %BeqzBeqv  to , size of [nBeqz, nBeqv]        
    end
    %BeqvBeqv
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqvPertBeqv evaluated in x+pert
        [dSf_dBeqv_PertBeqv, dSt_dBeqv_PertBeqv] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertBeqv
        %2nd Derivatives of Sbr w.r.t. BeqvBeqv
        d2Sf_dBeqv2(:, k) = (dSf_dBeqv_PertBeqv - dSf_dBeqv).' * lam / pert;  %BeqvBeqv from, size of [nBeqv, nBeqv]
        d2St_dBeqv2(:, k) = (dSt_dBeqv_PertBeqv - dSt_dBeqv).' * lam / pert;  %BeqvBeqv  to , size of [nBeqv, nBeqv]
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf14 = sparse(d2Sf_dBeqvVa);
num_Hf24 = sparse(d2Sf_dBeqvVm);
num_Hf34 = sparse(d2Sf_dBeqvBeqz);

num_Hf41 = sparse(d2Sf_dVaBeqv);
num_Hf42 = sparse(d2Sf_dVmBeqv);
num_Hf43 = sparse(d2Sf_dBeqzBeqv);

num_Hf44 = sparse(d2Sf_dBeqv2);

num_Ht14 = sparse(d2St_dBeqvVa);
num_Ht24 = sparse(d2St_dBeqvVm);
num_Ht34 = sparse(d2St_dBeqvBeqz);

num_Ht41 = sparse(d2St_dVaBeqv);
num_Ht42 = sparse(d2St_dVmBeqv);
num_Ht43 = sparse(d2St_dBeqzBeqv);

num_Ht44 = sparse(d2St_dBeqv2);


