function [num_G16, num_G26, num_G56, num_G61, num_G62, num_G65, num_G66] = d2Sbus_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXBEQV2PERT   Computes 2nd derivatives of power injection w.r.t. BeqvVa, BeqvVm, BeqvBz, VaBeqv, VmBeqv, BeqzBeqv, BeqvBeqv (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = D2SBUS_DXBEQV2(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)
%   [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = D2SBUS_DXBEQV2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 7 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. BeqvVa, BeqvVm, BeqvBz, VaBeqv, VmBeqv, BeqzBeqv,
%   and BeqvBeqv
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = d2Sbus_dxBeqv2(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GBvVa = d/dVa   (dSbus_dBeqv.' * lam)
%           GBvVm = d/dVm   (dSbus_dBeqv.' * lam)
%           GBvBz = d/dBeqz (dSbus_dBeqv.' * lam)
%           GVaBv = d/dBeqv (dSbus_dVa.'   * lam)
%           GVmBv = d/dBeqv (dSbus_dVm.'   * lam)
%           GBzBv = d/dBeqv (dSbus_dBeqz.' * lam)
%           GBvBv = d/dBeqv (dSbus_dBeqv.' * lam)

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
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC size[nBeqv,1]
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

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
    [dSbus_dBeqv] = dSbus_dBeq(branch, V, 2, vcart);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the selector
    BeqzAux1(iBeqz) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
            
    %Selector of active Beqv 
    BeqvAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the selector
    BeqvAux1(iBeqv) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqvAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed)  
    d2Sbus_dBeqvVa   = sparse( zeros(nb   ,nBeqv) ); 
    d2Sbus_dBeqvVm   = sparse( zeros(nb   ,nBeqv) ); 
    d2Sbus_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    
    d2Sbus_dVaBeqv   = sparse( zeros(nBeqv,nb   ) ); 
    d2Sbus_dVmBeqv   = sparse( zeros(nBeqv,nb   ) ); 
    d2Sbus_dBeqzBeqv = sparse( zeros(nBeqv,nBeqz) );
    
    d2Sbus_dBeqv2    = sparse( zeros(nBeqv,nBeqv) );

    %BeqvVa
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dBeqv_PertVa] = dSbus_dBeq(branch, V1p, 2, vcart); %dSbus_dBeqvPertVa
        d2Sbus_dVaBeqv(:, k) = (dSbus_dBeqv_PertVa - dSbus_dBeqv).' * lam / pert; %BeqvVa (dSbus_dBeqvPertVa - dSbus_dBeqv) size of [nBeqv, nb] 
    end
    %BeqvVm
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dBeqv_PertVm] = dSbus_dBeq(branch, V2p, 2, vcart); %dSbus_dBeqvPertVm
        d2Sbus_dVmBeqv(:, k) = (dSbus_dBeqv_PertVm - dSbus_dBeqv).' * lam / pert; %BeqvVm (dSbus_dBeqvPertVm - dSbus_dBeqv) size of [nBeqv, nb] 
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
        %dSbus_dBeqvPertBeqz evaluated in x+pert
        [dSbus_dBeqv_PertBeqz] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertBeqz
        %2nd Derivatives of Sbus w.r.t. BeqvBeqz
        d2Sbus_dBeqzBeqv(:, k) = (dSbus_dBeqv_PertBeqz - dSbus_dBeqv).' * lam / pert;  %BeqvBeqz (dSbus_dBeqvPertBeqz - dSbus_dBeqv) size of [nBeqv, nBeqz] 
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
        [dSbus_dV1_PertBeqv, dSbus_dV2_PertBeqv] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. BeqvVx
        d2Sbus_dBeqvVa(:, k) = (dSbus_dV1_PertBeqv - dSbus_dV1).' * lam / pert;  %VaBeqv (dSbus_dVaPertBeqv - dSbus_dVa) %AAB- Final d2Sbus_dBeqVa has a size of [nb, nBeqx] 
        d2Sbus_dBeqvVm(:, k) = (dSbus_dV2_PertBeqv - dSbus_dV2).' * lam / pert;  %VmBeqv (dSbus_dVmPertBeqv - dSbus_dVm) %AAB- Final d2Sbus_dBeqVm has a size of [nb, nBeqx] 
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
        %dSbus_dBeqvPertBeqv evaluated in x+pert
        [dSbus_dBeqz_PertBeqv] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertBeqzv
        %2nd Derivatives of Sbus w.r.t. BeqzBeqv
        d2Sbus_dBeqvBeqz(:, k) = (dSbus_dBeqz_PertBeqv - dSbus_dBeqz).' * lam / pert;  %BeqzBeqv (dSbus_dBeqzPertBeqv - dSbus_dBeqz) size of [nBeqz, nBeqv] 
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
        %dSbus_dBeqvPertBeqv evaluated in x+pert
        [dSbus_dBeqv_PertBeqv] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertBeqv
        %2nd Derivatives of Sbus w.r.t. Beqv2   
        d2Sbus_dBeqv2(:, k) = (dSbus_dBeqv_PertBeqv - dSbus_dBeqv).' * lam / pert;  %BeqvBeqv (dSbus_dBeqvPertBeqv - dSbus_dBeqv) %AAB- Final d2Sbus_dBeq2 has a size of [nBeqv , nBeqv] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G16 = sparse(d2Sbus_dBeqvVa);
num_G26 = sparse(d2Sbus_dBeqvVm);
num_G56 = sparse(d2Sbus_dBeqvBeqz);
num_G61 = sparse(d2Sbus_dVaBeqv);
num_G62 = sparse(d2Sbus_dVmBeqv);
num_G65 = sparse(d2Sbus_dBeqzBeqv);
num_G66 = sparse(d2Sbus_dBeqv2);
