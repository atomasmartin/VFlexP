function [num_G18, num_G28, num_G58, num_G68, num_G78, num_G81, num_G82, num_G85, num_G86, num_G87, num_G88] = d2Sbus_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SBUS_DXQTMA2PERT   Computes 2nd derivatives of power injection w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma (Finite Differences Method).
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = D2SBUS_DXQTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, VCART)

%   [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = D2SBUS_DXQTMA2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 11 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, 
%   GBeqzqtma, GBeqvqtma, GShqtma and Gqtmaqtma
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = ...
%                                                                           d2Sbus_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
%
%       Here the output matrices correspond to:
%           GqtmaVa = d/dVa   (dSbus_dqtma.' * lam)
%           GqtmaVm = d/dVm   (dSbus_dqtma.' * lam)
%           GqtmaBz = d/dBeqz (dSbus_dqtma.' * lam)
%           GqtmaBv = d/dBeqv (dSbus_dqtma.' * lam)
%           GqtmaSh = d/dsh   (dSbus_dqtma.' * lam)
%           GVaqtma = d/dqtma (dSbus_dVa.'   * lam)
%           GVmqtma = d/dqtma (dSbus_dVm.'   * lam)
%           GBzqtma = d/dqtma (dSbus_dBeqz.' * lam)
%           GBvqtma = d/dqtma (dSbus_dBeqv.' * lam)
%           GShqtma = d/dqtma (dSbus_dsh.'   * lam)
%           Gmaqtma = d/dqtma (dSbus_dqtma.' * lam)
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

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxqtma2Pert: Derivatives of Power balance equations w.r.t ma/tap Finite Differences Method in cartasian have not been coded yet')    

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
    diagShSel = sparse( diag(ShSel));         %AAB- Diagonal of the selector for derivative w.r.t. Theta_shift, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) );     %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Selector of active Qt ma/tap 
    QtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    QtmaSel(iQtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Qt
    diagQtmaSel = sparse( diag(QtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Qtma, size [nl,nl]
    diagYsQtma = sparse( diag(QtmaSel.*Ys) ); %AAB- Qtma selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)  
    d2Sbus_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Sbus_dQtmaVm   = sparse( zeros(nb,   nQtma) ); 
    d2Sbus_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Sbus_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Sbus_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) );
    
    d2Sbus_dVaQtma   = sparse( zeros(nQtma, nb   ) ); 
    d2Sbus_dVmQtma   = sparse( zeros(nQtma, nb   ) ); 
    d2Sbus_dBeqzQtma = sparse( zeros(nQtma, nBeqz) );
    d2Sbus_dBeqvQtma = sparse( zeros(nQtma, nBeqv) );
    d2Sbus_dPfshQtma = sparse( zeros(nQtma, nPfsh) ); 
    
    d2Sbus_dQtma2    = sparse( zeros(nQtma,nQtma) );
    
    %QtmaVa num_G81
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSbus_dQtma_PertVa] = dSbus_dma(branch, V1p, 2, vcart); %dSbus_dQtmaPertVa
        d2Sbus_dVaQtma(:, k) = (dSbus_dQtma_PertVa - dSbus_dQtma).' * lam / pert; %QtmaVa (dSbus_dQtmaPertVa - dSbus_dQtma) size of [nQtma, nb] 
    end
    %QtmaVm num_G82
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSbus_dQtma_PertVm] = dSbus_dma(branch, V2p, 2, vcart); %dSbus_dQtmaPertVm
        d2Sbus_dVmQtma(:, k) = (dSbus_dQtma_PertVm - dSbus_dQtma).' * lam / pert; %QtmaVm (dSbus_dQtmaPertVm - dSbus_dQtma) size of [nQtma, nb] 
    end
    %QtmaBeqz num_G85
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertBeqz evaluated in x+pert
        [dSbus_dQtma_PertBeqz] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertBeqz
        %2nd Derivatives of Sbus w.r.t. QtmaBeqz
        d2Sbus_dBeqzQtma(:, k) = (dSbus_dQtma_PertBeqz - dSbus_dQtma).' * lam / pert;  %QtmaBeqz (dSbus_dQtmaPertBeqz - dSbus_dQtma) size of [nQtma, nBeqz] 
    end
    %QtmaBeqv num_G86
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertBeqv evaluated in x+pert
        [dSbus_dQtma_PertBeqv] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertBeqv
        %2nd Derivatives of Sbus w.r.t. QtmaBeqv
        d2Sbus_dBeqvQtma(:, k) = (dSbus_dQtma_PertBeqv - dSbus_dQtma).' * lam / pert;  %QtmaBeqv (dSbus_dQtmaPertBeqv - dSbus_dQtma) size of [nQtma, nBeqv] 
    end
    %QtmaPfsh num_G87
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertPfsh evaluated in x+pert
        [dSbus_dQtma_PertPfsh] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertPfsh
        %2nd Derivatives of Sbus w.r.t. QtmaPfsh
        d2Sbus_dPfshQtma(:, k) = (dSbus_dQtma_PertPfsh - dSbus_dQtma).' * lam / pert;  %QtmaPfsh (dSbus_dQtmaPertPfsh - dSbus_dQtma) size of [nQtma, nPfsh] 
    end    
    
    %VxQtma num_G18 num_G28
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dVaPertQtma evaluated in x+pert
        [dSbus_dV1_PertQtma, dSbus_dV2_PertQtma] = dSbus_dV(Ybus_Pert, V, vcart);
        %2nd Derivatives of Sbus w.r.t. QtmaVx
        d2Sbus_dQtmaVa(:, k) = (dSbus_dV1_PertQtma - dSbus_dV1).' * lam / pert;  %VaQtma (dSbus_dVaPertQtma - dSbus_dVa) %AAB- Final d2Sbus_dQtmaVa has a size of [nb, nQtma] 
        d2Sbus_dQtmaVm(:, k) = (dSbus_dV2_PertQtma - dSbus_dV2).' * lam / pert;  %VmQtma (dSbus_dVmPertQtma - dSbus_dVm) %AAB- Final d2Sbus_dQtmaVm has a size of [nb, nQtma] 
    end
    %BeqzQtma num_G58
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqzPertQtma evaluated in x+pert
        [dSbus_dBeqz_PertQtma] = dSbus_dBeq(branch_Pert, V, 1, vcart); %dSbus_dBeqzPertQtma
        %2nd Derivatives of Sbus w.r.t. BeqzQtma
        d2Sbus_dQtmaBeqz(:, k) = (dSbus_dBeqz_PertQtma - dSbus_dBeqz).' * lam / pert;  %BeqzQtma (dSbus_dBeqzPertQtma - dSbus_dBeqz) size of [nBeqz, nQtma] 
    end
    %BeqvQtma num_G68
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dBeqvPertQtma evaluated in x+pert
        [dSbus_dBeqv_PertQtma] = dSbus_dBeq(branch_Pert, V, 2, vcart); %dSbus_dBeqvPertQtma
        %2nd Derivatives of Sbus w.r.t. BeqvQtma
        d2Sbus_dQtmaBeqv(:, k) = (dSbus_dBeqv_PertQtma - dSbus_dBeqv).' * lam / pert;  %BeqvQtma (dSbus_dBeqvPertQtma - dSbus_dBeqv) size of [nBeqv, nQtma] 
    end
    %PfshQtma num_G78
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dPfshPertQtma evaluated in x+pert
        [dSbus_dPfsh_PertQtma] = dSbus_dsh(branch_Pert, V, 1, vcart); %dSbus_dPfshPertQtma
        %2nd Derivatives of Sbus w.r.t. PfshQtma   
        d2Sbus_dQtmaPfsh(:, k) = (dSbus_dPfsh_PertQtma - dSbus_dPfsh).' * lam / pert;  %PfshQtma (dSbus_dPfshPertQtma - dSbus_dPfsh) size of [nPfsh , nQtma] 
    end
    %Qtma2 num_G88
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbus_dQtmaPertQtma evaluated in x+pert
        [dSbus_dQtma_PertQtma] = dSbus_dma(branch_Pert, V, 2, vcart); %dSbus_dQtmaPertQtma
        %2nd Derivatives of Sbus w.r.t. Qtma2   
        d2Sbus_dQtma2(:, k) = (dSbus_dQtma_PertQtma - dSbus_dQtma).' * lam / pert;  %Qtma2 (dSbus_dQtmaPertQtma - dSbus_dQtma) size of [nQtma , nQtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_G18 = sparse(d2Sbus_dQtmaVa);
num_G28 = sparse(d2Sbus_dQtmaVm);
num_G58 = sparse(d2Sbus_dQtmaBeqz);
num_G68 = sparse(d2Sbus_dQtmaBeqv);
num_G78 = sparse(d2Sbus_dQtmaPfsh);

num_G81 = sparse(d2Sbus_dVaQtma);
num_G82 = sparse(d2Sbus_dVmQtma);
num_G85 = sparse(d2Sbus_dBeqzQtma);
num_G86 = sparse(d2Sbus_dBeqvQtma);
num_G87 = sparse(d2Sbus_dPfshQtma);

num_G88 = sparse(d2Sbus_dQtma2);

