function [num_G19, num_G29, num_G59, num_G69, num_G79, num_G89, num_G91, num_G92, num_G95, num_G96, num_G97, num_G98, num_G99] = d2Sbus_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXVTMA2PERT   Computes 2nd derivatives of power injection w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma, qtmaqtma (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GvtmaVa, GvtmaVm, GvtmaBeqz, GvtmaBeqv, GvtmaSh, Gvtmaqtma, GVavtma, GVmvtma, GBeqzvtma, GBeqvvtma, GShvtma, Gqtmavtma, Gqtmaqtma] = D2SBUS_DXVTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)

%   [GvtmaVa, GvtmaVm, GvtmaBeqz, GvtmaBeqv, GvtmaSh, Gvtmaqtma, GVavtma, GVmvtma, GBeqzvtma, GBeqvvtma, GShvtma, Gqtmavtma, Gqtmaqtma] = D2SBUS_DXQTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 13 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, 
%   Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma and qtmaqtma.
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GvtmaVa, GvtmaVm, GvtmaBeqz, GvtmaBeqv, GvtmaSh, Gvtmaqtma, GVavtma,...
%                 GVmvtma, GBeqzvtma, GBeqvvtma, GShvtma, Gqtmavtma, Gqtmaqtma] = ...
%                                                             d2Sbus_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GvtmaVa   = d/dVa   (dSbus_dvtma.' * lam)
%           GvtmaVm   = d/dVm   (dSbus_dvtma.' * lam)
%           GvtmaBz   = d/dBeqz (dSbus_dvtma.' * lam)
%           GvtmaBv   = d/dBeqv (dSbus_dvtma.' * lam)
%           GvtmaSh   = d/dsh   (dSbus_dvtma.' * lam)
%           Gvtmavtma = d/dvtma (dSbus_dvtma.' * lam)
%           GVavtma   = d/dvtma (dSbus_dVa.'   * lam)
%           GVmvtma   = d/dvtma (dSbus_dVm.'   * lam)
%           GBzvtma   = d/dvtma (dSbus_dBeqz.' * lam)
%           GBvvtma   = d/dvtma (dSbus_dBeqv.' * lam)
%           GShvtma   = d/dvtma (dSbus_dsh.'   * lam)
%           Gqtmavtma = d/dvtma (dSbus_dqtma.' * lam)
%           Gvtmavtma = d/dvtma (dSbus_dvtma.' * lam)
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
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxvtma2Pert: Derivatives of Power balance equations w.r.t ma/tap Finite Differences Method in cartasian have not been coded yet')    

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
        
    %Dimensionalize (Allocate for computational speed)  
    d2Sbus_dVtmaVa   = sparse( zeros(nb,   nVtma) ); 
    d2Sbus_dVtmaVm   = sparse( zeros(nb,   nVtma) ); 
    d2Sbus_dVtmaBeqz = sparse( zeros(nBeqz,nVtma) );
    d2Sbus_dVtmaBeqv = sparse( zeros(nBeqv,nVtma) );
    d2Sbus_dVtmaPfsh = sparse( zeros(nPfsh,nVtma) );
    d2Sbus_dVtmaQtma = sparse( zeros(nQtma,nVtma) );
    
    d2Sbus_dVaVtma   = sparse( zeros(nVtma,nb   ) ); 
    d2Sbus_dVmVtma   = sparse( zeros(nVtma,nb   ) ); 
    d2Sbus_dBeqzVtma = sparse( zeros(nVtma,nBeqz) );
    d2Sbus_dBeqvVtma = sparse( zeros(nVtma,nBeqv) );
    d2Sbus_dPfshVtma = sparse( zeros(nVtma,nPfsh) );
    d2Sbus_dQtmaVtma = sparse( zeros(nVtma,nQtma) );    
    
    d2Sbus_dVtma2    = sparse( zeros(nVtma,nVtma) );
    
    %VtmaVa num_G91
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dVtma_PertVa] = dSbus_dma(branch, V1p, 4, vcart); %dSbus_dVtmaPertVa
        d2Sbus_dVaVtma(:, k) = (dSbus_dVtma_PertVa - dSbus_dVtma).' * lam / pert; %VtmaVa (dSbus_dVtmaPertVa - dSbus_dVtma) size of [nVtma, nb] 
    end
    %VtmaVm num_G92
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dVtma_PertVm] = dSbus_dma(branch, V2p, 4, vcart); %dSbus_dVtmaPertVm
        d2Sbus_dVmVtma(:, k) = (dSbus_dVtma_PertVm - dSbus_dVtma).' * lam / pert; %VtmaVm (dSbus_dVtmaPertVm - dSbus_dVtma) size of [nVtma, nb] 
    end
    %VtmaBeqz num_G95
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVtmaPertBeqz evaluated in x+pert
        [dSbus_dVtma_PertBeqz] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertBeqz
        %2nd Derivatives of Sbus w.r.t. VtmaBeqz
        d2Sbus_dBeqzVtma(:, k) = (dSbus_dVtma_PertBeqz - dSbus_dVtma).' * lam / pert;  %VtmaBeqz (dSbus_dVtmaPertBeqz - dSbus_dVtma) size of [nVtma, nBeqz] 
    end
    %VtmaBeqv num_G96
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVtmaPertBeqv evaluated in x+pert
        [dSbus_dVtma_PertBeqv] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertBeqv
        %2nd Derivatives of Sbus w.r.t. VtmaBeqv
        d2Sbus_dBeqvVtma(:, k) = (dSbus_dVtma_PertBeqv - dSbus_dVtma).' * lam / pert;  %VtmaBeqv (dSbus_dVtmaPertBeqv - dSbus_dVtma) size of [nVtma, nBeqv] 
    end
    %VtmaPfsh num_G97
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVtmaPertPfsh evaluated in x+pert
        [dSbus_dVtma_PertPfsh] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertPfsh
        %2nd Derivatives of Sbus w.r.t. VtmaPfsh
        d2Sbus_dPfshVtma(:, k) = (dSbus_dVtma_PertPfsh - dSbus_dVtma).' * lam / pert;  %VtmaPfsh (dSbus_dVtmaPertPfsh - dSbus_dVtma) size of [nVtma, nPfsh] 
    end    
    %VtmaQtma num_G98
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVtmaPertPfsh evaluated in x+pert
        [dSbus_dVtma_PertQtma] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertQtma
        %2nd Derivatives of Sbus w.r.t. VtmaQtma
        d2Sbus_dQtmaVtma(:, k) = (dSbus_dVtma_PertQtma - dSbus_dVtma).' * lam / pert;  %VtmaQtma (dSbus_dVtmaPertQtma - dSbus_dVtma) size of [nVtma, nQtma] 
    end 
    
    %VxVtma num_G19 num_G29
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVaPertVtma evaluated in x+pert
        [dSbus_dV1_PertVtma, dSbus_dV2_PertVtma] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. VtmaVx
        d2Sbus_dVtmaVa(:, k) = (dSbus_dV1_PertVtma - dSbus_dV1).' * lam / pert;  %VaVtma (dSbus_dVaPertVtma - dSbus_dVa) %AAB- Final d2Sbus_dVtmaVa has a size of [nb, nVtma] 
        d2Sbus_dVtmaVm(:, k) = (dSbus_dV2_PertVtma - dSbus_dV2).' * lam / pert;  %VmVtma (dSbus_dVmPertVtma - dSbus_dVm) %AAB- Final d2Sbus_dVtmaVm has a size of [nb, nVtma] 
    end
    %BeqzVtma num_G59
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqzPertVtma evaluated in x+pert
        [dSbus_dBeqz_PertVtma] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertVtma
        %2nd Derivatives of Sbus w.r.t. BeqzVtma
        d2Sbus_dVtmaBeqz(:, k) = (dSbus_dBeqz_PertVtma - dSbus_dBeqz).' * lam / pert;  %BeqzVtma (dSbus_dBeqzPertVtma - dSbus_dBeqz) size of [nBeqz, nVtma] 
    end
    %BeqvVtma num_G69
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqvPertVtma evaluated in x+pert
        [dSbus_dBeqv_PertVtma] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertVtma
        %2nd Derivatives of Sbus w.r.t. BeqvVtma
        d2Sbus_dVtmaBeqv(:, k) = (dSbus_dBeqv_PertVtma - dSbus_dBeqv).' * lam / pert;  %BeqvVtma (dSbus_dBeqvPertVtma - dSbus_dBeqv) size of [nBeqv, nVtma] 
    end
    %PfshVtma num_G79
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertVtma evaluated in x+pert
        [dSbus_dPfsh_PertVtma] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertVtma
        %2nd Derivatives of Sbus w.r.t. PfshVtma   
        d2Sbus_dVtmaPfsh(:, k) = (dSbus_dPfsh_PertVtma - dSbus_dPfsh).' * lam / pert;  %PfshVtma (dSbus_dPfshPertVtma - dSbus_dPfsh) size of [nPfsh , nVtma] 
    end
    %QtmaVtma num_G89
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertVtma evaluated in x+pert
        [dSbus_dQtma_PertVtma] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertVtma
        %2nd Derivatives of Sbus w.r.t. QtmaVtma   
        d2Sbus_dVtmaQtma(:, k) = (dSbus_dQtma_PertVtma - dSbus_dQtma).' * lam / pert;  %QtmaVtma2 (dSbus_dQtmaPertVtma - dSbus_dQtma) size of [nQtma , nVtma] 
    end
    %Vtma2 num_G99
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertVtma evaluated in x+pert
        [dSbus_dVtma_PertVtma] = dSbus_dma(branch_Pert, V, 4, vcart); %dSbus_dVtmaPertVtma
        %2nd Derivatives of Sbus w.r.t. Vtma2   
        d2Sbus_dVtma2(:, k) = (dSbus_dVtma_PertVtma - dSbus_dVtma).' * lam / pert;  %Vtma2 (dSbus_dVtmaPertVtma - dSbus_dVtma) size of [nVtma , nVtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G19 = sparse(d2Sbus_dVtmaVa);
num_G29 = sparse(d2Sbus_dVtmaVm);
num_G59 = sparse(d2Sbus_dVtmaBeqz);
num_G69 = sparse(d2Sbus_dVtmaBeqv);
num_G79 = sparse(d2Sbus_dVtmaPfsh);
num_G89 = sparse(d2Sbus_dVtmaQtma);

num_G91 = sparse(d2Sbus_dVaVtma);
num_G92 = sparse(d2Sbus_dVmVtma);
num_G95 = sparse(d2Sbus_dBeqzVtma);
num_G96 = sparse(d2Sbus_dBeqvVtma);
num_G97 = sparse(d2Sbus_dPfshVtma);
num_G98 = sparse(d2Sbus_dQtmaVtma);

num_G99 = sparse(d2Sbus_dVtma2);



