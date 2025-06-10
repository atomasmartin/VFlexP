function [num_G110, num_G210, num_G510, num_G610, num_G710, num_G810, num_G910, num_G101, num_G102, num_G105, num_G106, num_G107, num_G108, num_G109, num_G1010] = d2Sbus_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXSHDP2PERT Computes 2nd derivatives of power injection w.r.t. shdpVa, shdpVm, shdpBeqz, shdpBeqv, shdpSh, shdpqtma, shdpvtma, Vashdp, Vmshdp, Beqzshdp, Beqvshdp, Shshdp, qtmashdp, vtmashdp, shdpshdp, (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GshdpVa, GshdpVm, GshdpBeqz, GshdpBeqv, GshdpSh, Gshdpqtma, Gshdpvtma, GVashdp, GVmshdp, GBeqzshdp, GBeqvshdp, GShshdp, Gqtmashdp, Gvtmashdp, Gqtmashdp] = D2SBUS_DXSHDP2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)

%   [GshdpVa, GshdpVm, GshdpBeqz, GshdpBeqv, GshdpSh, Gshdpqtma, Gshdpvtma, GVashdp, GVmshdp, GBeqzshdp, GBeqvshdp, GShshdp, Gqtmashdp, Gvtmashdp, Gqtmashdp] = D2SBUS_DXSHDP2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 15 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. shdpVa, shdpVm, shdpBeqz, shdpBeqv, shdpSh, shdpqtma, shdpvtma, 
%   Vashdp, Vmshdp, Beqzshdp, Beqvshdp, Shshdp, qtmashdp, vtmashdp, and shdpshdp.
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GshdpVa, GshdpVm, GshdpBeqz, GshdpBeqv, GshdpSh, Gshdpqtma, Gshdpvtma,...
%                  GVashdp, GVmshdp, GBeqzshdp, GBeqvshdp, GShshdp, Gqtmashdp, Gvtmashdp, Gshdpshdp] = ...
%                                                             d2Sbus_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GshdpVa   = d/dVa   (dSbus_dshdp.' * lam)
%           GshdpVm   = d/dVm   (dSbus_dshdp.' * lam)
%           GshdpBz   = d/dBeqz (dSbus_dshdp.' * lam)
%           GshdpBv   = d/dBeqv (dSbus_dshdp.' * lam)
%           GshdpSh   = d/dsh   (dSbus_dshdp.' * lam)
%           Gshdpqtma = d/dqtma (dSbus_dshdp.' * lam)
%           Gshdpvtma = d/dvtma (dSbus_dshdp.' * lam)
%           GVashdp   = d/dshdp (dSbus_dVa.'   * lam)
%           GVmshdp   = d/dshdp (dSbus_dVm.'   * lam)
%           GBzshdp   = d/dshdp (dSbus_dBeqz.' * lam)
%           GBvshdp   = d/dshdp (dSbus_dBeqv.' * lam)
%           GShshdp   = d/dshdp (dSbus_dsh.'   * lam)
%           Gqtmashdp = d/dshdp (dSbus_dqtma.' * lam)
%           Gvtmashdp = d/dshdp (dSbus_dvtma.' * lam)
%           Gshdpshdp = d/dshdp (dSbus_dshdp.' * lam)
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
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxshdp2Pert: Derivatives of Power balance equations w.r.t Theta_dp Finite Differences Method in cartasian have not been coded yet')    

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
    [dSbus_dQtma] = dSbus_dma(branch, V, 2, vcart); 
    [dSbus_dVtma] = dSbus_dma(branch, V, 4, vcart); 
    [dSbus_dPfdp] = dSbus_dsh(branch, V, 3, vcart);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );         %AAB- Vector of zeros for the selector
    BeqzAux1(iBeqz) = 1;                      %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) );   %AAB- Beq Selector [nl,nBeqx]

    %Selector of active Beqv 
    BeqvAux1 = sparse( zeros(nl,1) );         %AAB- Vector of zeros for the selector
    BeqvAux1(iBeqv) = 1;                      %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) );   %AAB- Beq Selector [nl,nBeqx]

    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );            %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                         %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));         %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) );     %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Selector of active Qt ma/tap 
    QtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    QtmaSel(iQtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Qt
    diagQtmaSel = sparse( diag(QtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Qtma, size [nl,nl]
    diagYsQtma = sparse( diag(QtmaSel.*Ys) ); %AAB- Qtma selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Selector of active Vt ma/tap 
    VtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    VtmaSel(iVtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Vt
    diagVtmaSel = sparse( diag(VtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Vtma, size [nl,nl]
    diagYsVtma = sparse( diag(VtmaSel.*Ys) ); %AAB- Vtma selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Selector of active Theta_dp 
    ShdpSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    ShdpSel(iPfdp) = 1;                       %AAB- Fill the selector with 1 where Theta_dp is active controlling Pf
    diagShdpSel = sparse( diag(ShdpSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShdpShdp, size [nl,nl]
    diagYsshdp = sparse( diag(ShdpSel.*Ys) ); %AAB- Theta_dp selector multilied by the series addmitance Ys, size [nl,nl]
    
    
    %Dimensionalize (Allocate for computational speed)  
    d2Sbus_dPfdpVa   = sparse( zeros(nb,   nPfdp) ); 
    d2Sbus_dPfdpVm   = sparse( zeros(nb,   nPfdp) ); 
    d2Sbus_dPfdpBeqz = sparse( zeros(nBeqz,nPfdp) );
    d2Sbus_dPfdpBeqv = sparse( zeros(nBeqv,nPfdp) );
    d2Sbus_dPfdpPfsh = sparse( zeros(nPfsh,nPfdp) );
    d2Sbus_dPfdpQtma = sparse( zeros(nQtma,nPfdp) );
    d2Sbus_dPfdpVtma = sparse( zeros(nVtma,nPfdp) );
    
    d2Sbus_dVaPfdp   = sparse( zeros(nPfdp,nb   ) ); 
    d2Sbus_dVmPfdp   = sparse( zeros(nPfdp,nb   ) ); 
    d2Sbus_dBeqzPfdp = sparse( zeros(nPfdp,nBeqz) );
    d2Sbus_dBeqvPfdp = sparse( zeros(nPfdp,nBeqv) );
    d2Sbus_dPfshPfdp = sparse( zeros(nPfdp,nPfsh) );
    d2Sbus_dQtmaPfdp = sparse( zeros(nPfdp,nQtma) );    
    d2Sbus_dVtmaPfdp = sparse( zeros(nPfdp,nVtma) ); 
    
    d2Sbus_dPfdp2    = sparse( zeros(nPfdp,nPfdp) );
    
    %PfdpVa num_G101
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dPfdp_PertVa] = dSbus_dsh(branch, V1p, 3, vcart); %dSbus_dPfdpPertVa
        d2Sbus_dVaPfdp(:, k) = (dSbus_dPfdp_PertVa - dSbus_dPfdp).' * lam / pert; %PfdpVa (dSbus_dPfdpPertVa - dSbus_dPfdp) size of [nPfdp, nb] 
    end
    %PfdpVm num_G102
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dPfdp_PertVm] = dSbus_dsh(branch, V2p, 3, vcart); %dSbus_dPfdpPertVm
        d2Sbus_dVmPfdp(:, k) = (dSbus_dPfdp_PertVm - dSbus_dPfdp).' * lam / pert; %PfdpVm (dSbus_dPfdpPertVm - dSbus_dPfdp) size of [nPfdp, nb] 
    end
    %PfdpBeqz num_G105
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertBeqz evaluated in x+pert
        [dSbus_dPfdp_PertBeqz] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertBeqz
        %2nd Derivatives of Sbus w.r.t. PfdpBeqz
        d2Sbus_dBeqzPfdp(:, k) = (dSbus_dPfdp_PertBeqz - dSbus_dPfdp).' * lam / pert;  %PfdpBeqz (dSbus_dPfdpPertBeqz - dSbus_dPfdp) size of [nPfdp, nBeqz] 
    end
    %PfdpBeqv num_G106
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertBeqv evaluated in x+pert
        [dSbus_dPfdp_PertBeqv] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertBeqv
        %2nd Derivatives of Sbus w.r.t. PfdpBeqv
        d2Sbus_dBeqvPfdp(:, k) = (dSbus_dPfdp_PertBeqv - dSbus_dPfdp).' * lam / pert;  %PfdpBeqv (dSbus_dPfdpPertBeqv - dSbus_dPfdp) size of [nPfdp, nBeqv] 
    end
    %PfdpPfsh num_G107
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertPfsh evaluated in x+pert
        [dSbus_dPfdp_PertPfsh] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertPfsh
        %2nd Derivatives of Sbus w.r.t. PfdpPfsh
        d2Sbus_dPfshPfdp(:, k) = (dSbus_dPfdp_PertPfsh - dSbus_dPfdp).' * lam / pert;  %PfdpPfsh (dSbus_dPfdpPertPfsh - dSbus_dPfdp) size of [nPfdp, nPfsh] 
    end    
    %PfdpQtma num_G108
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertPfsh evaluated in x+pert
        [dSbus_dPfdp_PertQtma] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertQtma
        %2nd Derivatives of Sbus w.r.t. PfdpQtma
        d2Sbus_dQtmaPfdp(:, k) = (dSbus_dPfdp_PertQtma - dSbus_dPfdp).' * lam / pert;  %PfdpQtma (dSbus_dPfdpPertQtma - dSbus_dPfdp) size of [nPfdp, nQtma] 
    end 
    %PfdpVtma num_G109
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertPfsh evaluated in x+pert
        [dSbus_dPfdp_PertVtma] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertVtma
        %2nd Derivatives of Sbus w.r.t. PfdpVtma
        d2Sbus_dVtmaPfdp(:, k) = (dSbus_dPfdp_PertVtma - dSbus_dPfdp).' * lam / pert;  %PfdpVtma (dSbus_dPfdpPertVtma - dSbus_dPfdp) size of [nPfdp, nVtma] 
    end 
    
    %VxPfdp num_G110 num_G210
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVaPertShdp evaluated in x+pert
        [dSbus_dV1_PertPfdp, dSbus_dV2_PertPfdp] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. VtmaVx
        d2Sbus_dPfdpVa(:, k) = (dSbus_dV1_PertPfdp - dSbus_dV1).' * lam / pert;  %VaPfdp (dSbus_dVaPertPfdp - dSbus_dVa) %AAB- Final d2Sbus_dPfdpVa has a size of [nb, nPfdp] 
        d2Sbus_dPfdpVm(:, k) = (dSbus_dV2_PertPfdp - dSbus_dV2).' * lam / pert;  %VmPfdp (dSbus_dVmPertPfdp - dSbus_dVm) %AAB- Final d2Sbus_dPfdpVm has a size of [nb, nPfdp] 
    end
    %BeqzPfdp num_G510
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqzPertPfdp evaluated in x+pert
        [dSbus_dBeqz_PertPfdp] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertPfdp
        %2nd Derivatives of Sbus w.r.t. BeqzPfdp
        d2Sbus_dPfdpBeqz(:, k) = (dSbus_dBeqz_PertPfdp - dSbus_dBeqz).' * lam / pert;  %BeqzPfdp (dSbus_dBeqzPertPfdp - dSbus_dBeqz) size of [nBeqz, nPfdp] 
    end
    %BeqvVtma num_G610
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqvPertPfdp evaluated in x+pert
        [dSbus_dBeqv_PertPfdp] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertPfdp
        %2nd Derivatives of Sbus w.r.t. BeqvPfdp
        d2Sbus_dPfdpBeqv(:, k) = (dSbus_dBeqv_PertPfdp - dSbus_dBeqv).' * lam / pert;  %BeqvPfdp (dSbus_dBeqvPertPfdp - dSbus_dBeqv) size of [nBeqv, nPfdp] 
    end
    %PfshPfdp num_G710
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpSel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertVtma evaluated in x+pert
        [dSbus_dPfsh_PertPfdp] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertVtma
        %2nd Derivatives of Sbus w.r.t. PfshPfdp   
        d2Sbus_dPfdpPfsh(:, k) = (dSbus_dPfsh_PertPfdp - dSbus_dPfsh).' * lam / pert;  %PfshPfdp (dSbus_dPfshPertPfdp - dSbus_dPfsh) size of [nPfsh , nPfdp] 
    end
    %QtmaPfdp num_G810
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpSel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel);  
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertPfdp evaluated in x+pert
        [dSbus_dQtma_PertPfdp] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertPfdp
        %2nd Derivatives of Sbus w.r.t. QtmaPfdp   
        d2Sbus_dPfdpQtma(:, k) = (dSbus_dQtma_PertPfdp - dSbus_dQtma).' * lam / pert;  %QtmaPfdp2 (dSbus_dQtmaPertPfdp - dSbus_dQtma) size of [nQtma , nPfdp] 
    end
    %VtmaPfdp num_G910
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpSel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel);  
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVtmaPertPfdp evaluated in x+pert
        [dSbus_dVtma_PertPfdp] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertPfdp
        %2nd Derivatives of Sbus w.r.t. VtmaPfdp   
        d2Sbus_dPfdpVtma(:, k) = (dSbus_dVtma_PertPfdp - dSbus_dVtma).' * lam / pert;  %VtmaPfdp2 (dSbus_dVtmaPertPfdp - dSbus_dVtma) size of [nVtma , nPfdp] 
    end    
    %Pfdp2 num_G1010
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpSel representing only the active Pfdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfdpPertPfdp evaluated in x+pert
        [dSbus_dPfdp_PertPfdp] = dSbus_dsh(branch_Pert, V, 3, vcart); %dSbus_dPfdpPertPfdp
        %2nd Derivatives of Sbus w.r.t. Pfdp2   
        d2Sbus_dPfdp2(:, k) = (dSbus_dPfdp_PertPfdp - dSbus_dPfdp).' * lam / pert;  %Pfdp2 (dSbus_dPfdpPertPfdp - dSbus_dPfdp) size of [nPfdp , nPfdp] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G110 = sparse(d2Sbus_dPfdpVa);
num_G210 = sparse(d2Sbus_dPfdpVm);
num_G510 = sparse(d2Sbus_dPfdpBeqz);
num_G610 = sparse(d2Sbus_dPfdpBeqv);
num_G710 = sparse(d2Sbus_dPfdpPfsh);
num_G810 = sparse(d2Sbus_dPfdpQtma);
num_G910 = sparse(d2Sbus_dPfdpVtma);

num_G101 = sparse(d2Sbus_dVaPfdp);
num_G102 = sparse(d2Sbus_dVmPfdp);
num_G105 = sparse(d2Sbus_dBeqzPfdp);
num_G106 = sparse(d2Sbus_dBeqvPfdp);
num_G107 = sparse(d2Sbus_dPfshPfdp);
num_G108 = sparse(d2Sbus_dQtmaPfdp);
num_G109 = sparse(d2Sbus_dVtmaPfdp);

num_G1010 = sparse(d2Sbus_dPfdp2);
