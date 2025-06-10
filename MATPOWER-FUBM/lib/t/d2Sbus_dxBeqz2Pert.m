function [num_G15, num_G25, num_G51, num_G52, num_G55] = d2Sbus_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXBEQZ2PERT   Computes 2nd derivatives of power injection w.r.t. BeqzVa, BeqzVm, VaBeqz, VmBeqz, BeqzBeqz (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GBzVa, GBzVm, GVaBz, GVmBz, GBzBz] = D2SBUS_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT)
%   [GBzVa, GBzVm, GVaBz, GVmBz, GBzBz] = D2SBUS_DXBEQZ2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%
%   Returns 5 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. BeqzVa, BeqzVm, VaBeqz, VmBeqz, BeqzBeqz.
%
%
%   Examples:
%       [GBzVa, GBzVm, GVaBz, GVmBz, GBzBz] = d2Sbus_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GBzVa = d/dVa   (dSbus_dBeqz.' * lam)
%           GBzVm = d/dVm   (dSbus_dBeqz.' * lam)
%           GVaBz = d/dBeqz (dSbus_dVa.'   * lam)
%           GVmBz = d/dBeqz (dSbus_dVm.'   * lam)
%           GBzBz = d/dBeqz (dSbus_dBeq.'  * lam)
%
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
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
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
%[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control

if vcart
    error('d2Sbus_dxBeq2Pert: Derivatives of Power balance equations w.r.t Beq using Finite Differences Method in cartasian has not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Make Ybus, Yf, Yt
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
    %Sbus 1st Derivatives 
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart);
    [dSbus_dBeqz] = dSbus_dBeq(branch, V, 1, vcart);
    
    %Selector of active Beqz 
    BeqAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the seclector
    BeqAux1(iBeqz) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqsel = sparse( diag(BeqAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed)   
    d2Sbus_dVaBeqz = sparse( zeros(nBeqz,    nb) );
    d2Sbus_dVmBeqz = sparse( zeros(nBeqz,    nb) );
    d2Sbus_dBeqzVa = sparse( zeros(nb,    nBeqz) );
    d2Sbus_dBeqzVm = sparse( zeros(nb,    nBeqz) );
    d2Sbus_dBeqz2  = sparse( zeros(nBeqz, nBeqz) );
    
    %BeqzVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dBeqz_PertVa] = dSbus_dBeq(branch, V1p, 1, vcart); %dSbus_dBeqzPertVa
        d2Sbus_dVaBeqz(:, k) = (dSbus_dBeqz_PertVa - dSbus_dBeqz).' * lam / pert; %BeqzVa (dSbus_dBeqzPertVa - dSbus_dBeqz)
    end
    %BeqzVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dBeqz_PertVm] = dSbus_dBeq(branch, V2p, 1, vcart); %dSbus_dBeqzPertVm
        d2Sbus_dVmBeqz(:, k) = (dSbus_dBeqz_PertVm - dSbus_dBeqz).' * lam / pert; %BeqzVm (dSbus_dBeqzPertVm - dSbus_dBeqz)
    end
    %VxBeqz
    for k=1:nBeqz 
        PertSel=diagBeqsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVaPertBeqz evaluated in x+pert
        [dSbus_dV1_PertBeqz, dSbus_dV2_PertBeqz] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. BeqzVx
        d2Sbus_dBeqzVa(:, k) = (dSbus_dV1_PertBeqz - dSbus_dV1).' * lam / pert;  %VaBeqz (dSbus_dVaPertBeqz - dSbus_dVa) %AAB- Final d2Sbus_dBeqVa has a size of [nb, nBeqx] 
        d2Sbus_dBeqzVm(:, k) = (dSbus_dV2_PertBeqz - dSbus_dV2).' * lam / pert;  %VmBeqz (dSbus_dVmPertBeqz - dSbus_dVm) %AAB- Final d2Sbus_dBeqVm has a size of [nb, nBeqx] 
    end
    %BeqzBeqz
    for k=1:nBeqz 
        PertSel=diagBeqsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqzPertBeqz evaluated in x+pert
        [dSbus_dBeqz_PertBeqz] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertBeqz
        %2nd Derivatives of Sbus w.r.t. Beqz2   
        d2Sbus_dBeqz2(:, k) = (dSbus_dBeqz_PertBeqz - dSbus_dBeqz).' * lam / pert;  %BeqzBeqz (dSbus_dBeqzPertBeqz - dSbus_dBeqz) %AAB- Final d2Sbus_dBeq2 has a size of [n , n] 
        
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G15 = sparse(d2Sbus_dBeqzVa);
num_G25 = sparse(d2Sbus_dBeqzVm);
num_G51 = sparse(d2Sbus_dVaBeqz); %num_G15';
num_G52 = sparse(d2Sbus_dVmBeqz); %num_G25';
num_G55 = sparse(d2Sbus_dBeqz2);


