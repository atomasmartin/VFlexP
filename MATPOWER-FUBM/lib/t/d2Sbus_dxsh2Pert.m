function [num_G17, num_G27, num_G57, num_G67, num_G71, num_G72, num_G75, num_G76, num_G77] = d2Sbus_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXSH2PERT   Computes 2nd derivatives of power injection w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GShVa, GShVm, GShBz, GShBv, GVaSh, GVmSh, GBzSh, GBvSh, GShSh] = D2SBUS_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)

%   [GShVa, GShVm, GShBz, GShBv, GVaSh, GVmSh, GBzSh, GBvSh, GShSh] = D2SBUS_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 9 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. GShVa, GShVm, GShBz, GShBv, GVaSh, GVmSh, GBzSh, GBvSh,
%   and ShSh
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GShVa, GShVm, GShBz, GShBv, GVaSh, GVmSh, GBzSh, GBvSh, GShSh] = d2Sbus_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GShVa = d/dVa   (dSbus_dsh.'   * lam)
%           GShVm = d/dVm   (dSbus_dsh.'   * lam)
%           GShBz = d/dBeqz (dSbus_dsh.'   * lam)
%           GShBv = d/dBeqv (dSbus_dsh.'   * lam)
%           GVaSh = d/dsh   (dSbus_dVa.'   * lam)
%           GVmSh = d/dsh   (dSbus_dVm.'   * lam)
%           GBzSh = d/dsh   (dSbus_dBeqz.' * lam)
%           GBvSh = d/dsh   (dSbus_dBeqv.' * lam)
%           GShSh = d/dsh   (dSbus_dsh.'   * lam)

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

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxsh2Pert: Derivatives of Power balance equations w.r.t Theta_sh using Finite Differences Method in cartasian have not been coded yet')    

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
    [dSbus_dPfsh] = dSbus_dsh(branch, V, 1, vcart);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the selector
    BeqzAux1(iBeqz) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    
    %Selector of active Beqv 
    BeqvAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the selector
    BeqvAux1(iBeqv) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]

    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]

    %Dimensionalize (Allocate for computational speed)   
    d2Sbus_dPfshVa   = sparse( zeros(nb   ,nPfsh) ); 
    d2Sbus_dPfshVm   = sparse( zeros(nb   ,nPfsh) ); 
    d2Sbus_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2Sbus_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );
    
    d2Sbus_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2Sbus_dVmPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2Sbus_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2Sbus_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
    
    d2Sbus_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );
    
    %PfShVa num_G71
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dPfsh_PertVa] = dSbus_dsh(branch, V1p, 1, vcart); %dSbus_dPfshPertVa
        d2Sbus_dVaPfsh(:, k) = (dSbus_dPfsh_PertVa - dSbus_dPfsh).' * lam / pert; %shVa (dSbus_dPfshPertVa - dSbus_dPfsh) size of [nPfsh, nb] 
    end
    %PfshVm num_G72
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dPfsh_PertVm] = dSbus_dsh(branch, V2p, 1, vcart); %dSbus_dPfshPertVm
        d2Sbus_dVmPfsh(:, k) = (dSbus_dPfsh_PertVm - dSbus_dPfsh).' * lam / pert; %PfshVm (dSbus_dPfshPertVm - dSbus_dPfsh) size of [nPfsh, nb] 
    end
    %PfshBeqz num_G75
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertBeqz evaluated in x+pert
        [dSbus_dPfsh_PertBeqz] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertBeqz
        %2nd Derivatives of Sbus w.r.t. PfshBeqz
        d2Sbus_dBeqzPfsh(:, k) = (dSbus_dPfsh_PertBeqz - dSbus_dPfsh).' * lam / pert;  %PfshBeqz (dSbus_dPfshPertBeqz - dSbus_dPfsh) size of [nPfsh, nBeqz] 
    end
    %PfshBeqv num_G76
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertBeqv evaluated in x+pert
        [dSbus_dPfsh_PertBeqv] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertBeqv
        %2nd Derivatives of Sbus w.r.t. PfshBeqv
        d2Sbus_dBeqvPfsh(:, k) = (dSbus_dPfsh_PertBeqv - dSbus_dPfsh).' * lam / pert;  %PfshBeqv (dSbus_dPfshPertBeqv - dSbus_dPfsh) size of [nPfsh, nBeqv] 
    end
    %VxPfsh num_G17 num_G27
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVaPertPfsh evaluated in x+pert
        [dSbus_dV1_PertPfsh, dSbus_dV2_PertPfsh] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. PfshVx
        d2Sbus_dPfshVa(:, k) = (dSbus_dV1_PertPfsh - dSbus_dV1).' * lam / pert;  %VaPfsh (dSbus_dVaPertPfsh - dSbus_dVa) %AAB- Final d2Sbus_dPfshVa has a size of [nb, nPfsh] 
        d2Sbus_dPfshVm(:, k) = (dSbus_dV2_PertPfsh - dSbus_dV2).' * lam / pert;  %VmPfsh (dSbus_dVmPertPfsh - dSbus_dVm) %AAB- Final d2Sbus_dPfshVm has a size of [nb, nPfsh] 
    end
    %BeqzPfsh num_G57
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqzPertPfsh evaluated in x+pert
        [dSbus_dBeqz_PertPfsh] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertPfsh
        %2nd Derivatives of Sbus w.r.t. BeqzPfsh
        d2Sbus_dPfshBeqz(:, k) = (dSbus_dBeqz_PertPfsh - dSbus_dBeqz).' * lam / pert;  %BeqzPfsh (dSbus_dBeqzPertPfsh - dSbus_dBeqz) size of [nBeqz, nPfsh] 
    end
    %BeqvPfsh num_G67
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqvPertPfsh evaluated in x+pert
        [dSbus_dBeqv_PertPfsh] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertPfsh
        %2nd Derivatives of Sbus w.r.t. BeqvPfsh
        d2Sbus_dPfshBeqv(:, k) = (dSbus_dBeqv_PertPfsh - dSbus_dBeqv).' * lam / pert;  %BeqvPfsh (dSbus_dBeqvPertPfsh - dSbus_dBeqv) size of [nBeqv, nPfsh] 
    end
    %PfshPfsh num_G77
    for k=1:nPfsh
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShSel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertPfsh evaluated in x+pert
        [dSbus_dPfsh_PertPfsh] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertPfsh
        %2nd Derivatives of Sbus w.r.t. Pfsh2   
        d2Sbus_dPfsh2(:, k) = (dSbus_dPfsh_PertPfsh - dSbus_dPfsh).' * lam / pert;  %PfshPfsh (dSbus_dPfshPertPfsh - dSbus_dPfsh) %AAB- Final d2Sbus_dsh2 has a size of [nPfsh , nPfsh] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G17 = sparse(d2Sbus_dPfshVa);
num_G27 = sparse(d2Sbus_dPfshVm);
num_G57 = sparse(d2Sbus_dPfshBeqz);
num_G67 = sparse(d2Sbus_dPfshBeqv);

num_G71 = sparse(d2Sbus_dVaPfsh);
num_G72 = sparse(d2Sbus_dVmPfsh);
num_G75 = sparse(d2Sbus_dBeqzPfsh);
num_G76 = sparse(d2Sbus_dBeqvPfsh);

num_G77 = sparse(d2Sbus_dPfsh2);
