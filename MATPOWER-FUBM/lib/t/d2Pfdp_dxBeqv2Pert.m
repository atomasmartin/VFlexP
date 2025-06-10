function [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44] = d2Pfdp_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXBEQV2PERT  Computes 2nd derivatives of Droop Control w.r.t. BeqvVa, BeqvVm, BeqvBeqz, VaBeqv, VmBeqv, BeqzBeqv, BeqvBeqv (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv] = D2PFDP_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 7 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv] = D2PFDP_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfBvVa, HfBvVm, HfBvBz, HfVaBv, HfVmBv, HfBzBv, HfBvBv] = D2PFDP_DXBEQV2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%       Here the output matrices correspond to:
%           HfBvVa = d/dBeqv (dPfdp_dVa.'   * mu)
%           HfBvVm = d/dBeqv (dPfdp_dVm.'   * mu)
%           HfBvBz = d/dBeqv (dPfdp_dBeqz.' * mu)
%           HfVaBv = d/dVa   (dPfdp_dBeqv.' * mu)
%           HfVmBv = d/dVm   (dPfdp_dBeqv.' * mu)
%           HfBzBv = d/dBeqz (dPfdp_dBeqv.' * mu)
%           HfBvBv = d/dBeqv (dPfdp_dBeqv.' * mu)
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

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

%% Calculation of derivatives
if vcart
    error('d2Pfdp_dxBeqv2Pert: Derivatives of Flow Limit equations w.r.t Beq using Finite Differences Method in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
   
    %Make Ybus, Yf, Yt
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

    %Derivatives of Voltage Magnitude w.r.t. Voltage magnitude
    dVmf_dVm = sparse(zeros(nl,nb));                          % Initialize for speed [nl,nb]
    fdp = branch(iPfdp, F_BUS);                               % List of "from" buses with Voltage Droop Control [nPfdp, 1]
    Cfdp = sparse(1:nPfdp, fdp, ones(nPfdp, 1), nPfdp, nb);   % connection matrix for line & from buses with Voltage Droop Control [nPfdp, nb]
    dVmf_dVm(iPfdp,:)=Cfdp;        
    
    %Sbr 1st Derivatives 
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1, vcart);
    [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2, vcart);
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    dPfdp_dBeqv = -real(dSf_dBeqv);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
      
    %Selector of active Beqv 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqvAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed) 
    d2Pfdp_dBeqvVa   = sparse( zeros(nb   ,nBeqv) );
    d2Pfdp_dBeqvVm   = sparse( zeros(nb   ,nBeqv) );
    d2Pfdp_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    
    d2Pfdp_dVaBeqv   = sparse( zeros(nBeqv,   nb) );
    d2Pfdp_dVmBeqv   = sparse( zeros(nBeqv,   nb) ); 
    d2Pfdp_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );

    d2Pfdp_dBeqv2    = sparse( zeros(nBeqv,nBeqv) );
    
    %BeqvVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dBeqv_PertVa, dSt_dBeqv_PertVa] = dSbr_dBeq(branch, V1p, 2, vcart); %dSbr_dBeqvPertVa
        dPfdp_dBeqv_PertVa = -real(dSf_dBeqv_PertVa);
        d2Pfdp_dVaBeqv(:, k) = (dPfdp_dBeqv_PertVa - dPfdp_dBeqv).' * lam / pert; %BeqvVa From
    end
    %BeqzVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dBeqv_PertVm, dSt_dBeqv_PertVm] = dSbr_dBeq(branch, V2p, 2, vcart); %dSbr_dBeqvPertVm
        dPfdp_dBeqv_PertVm = -real(dSf_dBeqv_PertVm);
        d2Pfdp_dVmBeqv(:, k) = (dPfdp_dBeqv_PertVm - dPfdp_dBeqv).' * lam / pert; %BeqvVm From
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
        %dSf_dBeqvPertBeqz evaluated in x+pert
        [dSf_dBeqv_PertBeqz, dSt_dBeqv_PertBeqz] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertBeqz
        %dPfdp_dBeqvPertBeqz evaluated in x+pert
        dPfdp_dBeqv_PertBeqz = -real(dSf_dBeqv_PertBeqz);
        %2nd Derivatives of Pfdp w.r.t. BeqvBeqz
        d2Pfdp_dBeqzBeqv(:, k) = (dPfdp_dBeqv_PertBeqz - dPfdp_dBeqv).' * lam / pert; %BeqvBeqz From
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
        
        %dPfdp_dVxPertBeqv evaluated in x+pert
        dPfdp_dV1_PertBeqv = -real(dSf_dV1_PertBeqv);
        dPfdp_dV2_PertBeqv = -real(dSf_dV2_PertBeqv) + Kdp.*( dVmf_dVm );
                
        %2nd Derivatives of Sbus w.r.t. BeqvVx
        d2Pfdp_dBeqvVa(:, k) = (dPfdp_dV1_PertBeqv - dPfdp_dV1).' * lam / pert;  %VaBeqv from, size of [nb, nBeqv] 
        d2Pfdp_dBeqvVm(:, k) = (dPfdp_dV2_PertBeqv - dPfdp_dV2).' * lam / pert;  %VmBeqv from, size of [nb, nBeqv] 
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
        %dPfdp_dBeqzPertBeqv evaluated in x+pert
        dPfdp_dBeqz_PertBeqv = -real(dSf_dBeqz_PertBeqv);       
        
        %2nd Derivatives of Sbr w.r.t. BeqzBeqv
        d2Pfdp_dBeqvBeqz(:, k) = (dPfdp_dBeqz_PertBeqv - dPfdp_dBeqz).' * lam / pert;  %BeqzBeqv from, size of [nBeqz, nBeqv]
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
        %dPfdp_dBeqzPertBeqv evaluated in x+pert
        dPfdp_dBeqv_PertBeqv = -real(dSf_dBeqv_PertBeqv);  
        
        %2nd Derivatives of Sbr w.r.t. BeqvBeqv
        d2Pfdp_dBeqv2(:, k) = (dPfdp_dBeqv_PertBeqv - dPfdp_dBeqv).' * lam / pert;  %BeqvBeqv from, size of [nBeqv, nBeqv]
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf14 = sparse(d2Pfdp_dBeqvVa);
num_Hf24 = sparse(d2Pfdp_dBeqvVm);
num_Hf34 = sparse(d2Pfdp_dBeqvBeqz);

num_Hf41 = sparse(d2Pfdp_dVaBeqv);
num_Hf42 = sparse(d2Pfdp_dVmBeqv);
num_Hf43 = sparse(d2Pfdp_dBeqzBeqv);

num_Hf44 = sparse(d2Pfdp_dBeqv2);
