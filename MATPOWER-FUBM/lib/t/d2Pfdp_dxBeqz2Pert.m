function [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33] = d2Pfdp_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXBEQZ2PERT  Computes 2nd derivatives of Droop Control w.r.t. BeqzVa, BeqzVm, VaBeqz, VmBeqz, BeqzBeqz (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%
%   [HfBzVa, HfBzVm, HfVaBz, HfVmBz, HfBzBz] = D2PFDP_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 5 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HfBzVa, HfBzVm, HfVaBz, HfVmBz, HfBzBz] = D2PFDP_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%       [HfBzVa, HfBzVm, HfVaBz, HfVmBz, HfBzBz] = d2Pfdp_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%
%       Here the output matrices correspond to:
%           HfBzVa = d/dBeqz (dPfdp_dVa.'   * mu)
%           HfBzVm = d/dBeqz (dPfdp_dVm.'   * mu)
%           HfVaBz = d/dVa   (dPfdp_dBeqz.' * mu)
%           HfVmBz = d/dVm   (dPfdp_dBeqz.' * mu)
%           HfBzBz = d/dBeqz (dPfdp_dBeqz.' * mu)
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

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

if vcart
    error('d2Pfdp_dxBeqz2Pert: Derivatives of Flow Limit equations w.r.t Beq using Finite Differences Method in cartasian have not been coded yet')    

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
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed) 
    d2Pfdp_dBeqzVa = sparse( zeros(nb,   nBeqz) );
    d2Pfdp_dBeqzVm = sparse( zeros(nb,   nBeqz) );

    d2Pfdp_dVaBeqz = sparse( zeros(nBeqz,   nb) );
    d2Pfdp_dVmBeqz = sparse( zeros(nBeqz,   nb) ); 

    d2Pfdp_dBeqz2  = sparse( zeros(nBeqz,nBeqz) );
    
    %BeqzVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dBeqz_PertVa, dSt_dBeqz_PertVa] = dSbr_dBeq(branch, V1p, 1, vcart); %dSbr_dBeqzPertVa
        dPfdp_dBeqz_PertVa = -real(dSf_dBeqz_PertVa);
        d2Pfdp_dVaBeqz(:, k) = (dPfdp_dBeqz_PertVa - dPfdp_dBeqz).' * lam / pert; %BeqzVa From
    end
    %BeqzVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dBeqz_PertVm, dSt_dBeqz_PertVm] = dSbr_dBeq(branch, V2p, 1, vcart); %dSbr_dBeqzPertVm
        dPfdp_dBeqz_PertVm = -real(dSf_dBeqz_PertVm);
        d2Pfdp_dVmBeqz(:, k) = (dPfdp_dBeqz_PertVm - dPfdp_dBeqz).' * lam / pert; %BeqzVm From
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
        [dSf_dV1_PertBeqz, dSf_dV2_PertBeqz, dSt_dV1_PertBeqz, dSt_dV2_PertBeqz, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        
        %dPfdp_dVxPertBeqz evaluated in x+pert
        dPfdp_dV1_PertBeqz = -real(dSf_dV1_PertBeqz);
        dPfdp_dV2_PertBeqz = -real(dSf_dV2_PertBeqz) + Kdp.*( dVmf_dVm );
        
        %2nd Derivatives of Pfdp w.r.t. BeqzVx
        d2Pfdp_dBeqzVa(:, k) = (dPfdp_dV1_PertBeqz - dPfdp_dV1).' * lam / pert;  %VaBeqz from, size of [nb, nBeqz] 
        d2Pfdp_dBeqzVm(:, k) = (dPfdp_dV2_PertBeqz - dPfdp_dV2).' * lam / pert;  %VmBeqz from, size of [nb, nBeqz]     
    end
    %BeqzBeqz
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertBeqz evaluated in x+pert
        [dSf_dBeqz_PertBeqz, dSt_dBeqz_PertBeqz] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertBeqz
        %dPfdp_dBeqzPertBeqz evaluated in x+pert
        dPfdp_dBeqz_PertBeqz = -real(dSf_dBeqz_PertBeqz);
        
        %2nd Derivatives of Pfdp w.r.t. Beqz2   
        d2Pfdp_dBeqz2(:, k) = (dPfdp_dBeqz_PertBeqz - dPfdp_dBeqz).' * lam / pert;  %BeqzBeqz from, size of [nBeqz , nBeqz] 
    end

end
%Assigning the partial derivatives with their respective outputs. 
num_Hf13 = sparse(d2Pfdp_dBeqzVa);
num_Hf23 = sparse(d2Pfdp_dBeqzVm);

num_Hf31 = sparse(d2Pfdp_dVaBeqz);
num_Hf32 = sparse(d2Pfdp_dVmBeqz);

num_Hf33 = sparse(d2Pfdp_dBeqz2);
