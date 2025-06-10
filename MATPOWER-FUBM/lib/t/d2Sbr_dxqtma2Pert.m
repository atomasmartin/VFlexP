function [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66,...
          num_Ht16, num_Ht26, num_Ht36, num_Ht46, num_Ht56, num_Ht61, num_Ht62, num_Ht63, num_Ht64, num_Ht65, num_Ht66] = d2Sbr_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SF_DXQTMA2PERT  Computes 2nd derivatives of complex brch power flow "from" w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma,...
%    HtqtmaVa, HtqtmaVm, HtqtmaBz, HtqtmaBv, HtqtmaSh, HtVaqtma, HtVmqtma, HtBzqtma, HtBvqtma, HtShqtma, Htqtmaqtma] = D2SBR_DXQTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 22 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector MU with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma,...
%    HtqtmaVa, HtqtmaVm, HtqtmaBz, HtqtmaBv, HtqtmaSh, HtVaqtma, HtVmqtma, HtBzqtma, HtBvqtma, HtShqtma, Htqtmaqtma] = D2SBR_DXQTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma,...
%    HtqtmaVa, HtqtmaVm, HtqtmaBz, HtqtmaBv, HtqtmaSh, HtVaqtma, HtVmqtma, HtBzqtma, HtBvqtma, HtShqtma, Htqtmaqtma] = d2Sbr_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%
%       Here the output matrices correspond to:
%           HqtmaVa   = d/dqtma (dSf_dVa.'   * mu)
%           HqtmaVm   = d/dqtma (dSf_dVm.'   * mu)
%           HqtmaBz   = d/dqtma (dSf_dBeqz.' * mu)
%           HqtmaBv   = d/dqtma (dSf_dBeqv.' * mu)
%           HqtmaSh   = d/dqtma (dSf_dSh.'   * mu)
%           HVaqtma   = d/dVa   (dSf_dqtma.' * mu)
%           HVmqtma   = d/dVm   (dSf_dqtma.' * mu)
%           HBzqtma   = d/dBeqz (dSf_dqtma.' * mu)
%           HBvqtma   = d/dBeqv (dSf_dqtma.' * mu)
%           HShqtma   = d/dSh   (dSf_dqtma.' * mu)
%           Hqtmaqtma = d/dqtma (dSf_dqtma.' * mu)
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
%Changing Perturbation to Degrees for Theta_sh since inside makeYbus it gets changed to radians.
pertDeg = (pert*180)/pi;
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap

%% Calculation of derivatives
if vcart
    error('d2Sbr_dxma2Pert: Derivatives of Flow Limit equations w.r.t ma using Finite Differences Method in cartasian have not been coded yet')    

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
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, vcart);
    
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
        
    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
 
    %Selector of active Qt ma/tap 
    QtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    QtmaSel(iQtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Qt
    diagQtmaSel = sparse( diag(QtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Qtma, size [nl,nl]
    diagYsQtma = sparse( diag(QtmaSel.*Ys) ); %AAB- Qtma selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)  
    d2Sf_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Sf_dQtmaVm   = sparse( zeros(nb,   nQtma) );
    d2Sf_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Sf_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Sf_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) ); 

    d2Sf_dVaQtma   = sparse( zeros(nQtma,nb   ) ); 
    d2Sf_dVmQtma   = sparse( zeros(nQtma,nb   ) );
    d2Sf_dBeqzQtma = sparse( zeros(nQtma,nBeqz) );
    d2Sf_dBeqvQtma = sparse( zeros(nQtma,nBeqv) );
    d2Sf_dPfshQtma = sparse( zeros(nQtma,nPfsh) ); 
    
    d2Sf_dQtma2    = sparse( zeros(nQtma,nQtma) );

    d2St_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2St_dQtmaVm   = sparse( zeros(nb,   nQtma) );
    d2St_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2St_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2St_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) ); 

    d2St_dVaQtma   = sparse( zeros(nQtma,nb   ) ); 
    d2St_dVmQtma   = sparse( zeros(nQtma,nb   ) );
    d2St_dBeqzQtma = sparse( zeros(nQtma,nBeqz) );
    d2St_dBeqvQtma = sparse( zeros(nQtma,nBeqv) );
    d2St_dPfshQtma = sparse( zeros(nQtma,nPfsh) ); 
    
    d2St_dQtma2    = sparse( zeros(nQtma,nQtma) );    
    
    %QtmaVa num_G61
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dQtma_PertVa, dSt_dQtma_PertVa] = dSbr_dma(branch, V1p, 2, vcart); %dSbr_dQtmaPertVa
        d2Sf_dVaQtma(:, k) = (dSf_dQtma_PertVa - dSf_dQtma).' * lam / pert; %QtmaVa from, size of [nQtma, nb] 
        d2St_dVaQtma(:, k) = (dSt_dQtma_PertVa - dSt_dQtma).' * lam / pert; %QtmaVa  to , size of [nQtma, nb] 
    end
    %QtmaVm num_G62
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dQtma_PertVm, dSt_dQtma_PertVm] = dSbr_dma(branch, V2p, 2, vcart); %dSbr_dQtmaPertVm
        d2Sf_dVmQtma(:, k) = (dSf_dQtma_PertVm - dSf_dQtma).' * lam / pert; %QtmaVm from, size of [nQtma, nb] 
        d2St_dVmQtma(:, k) = (dSt_dQtma_PertVm - dSt_dQtma).' * lam / pert; %QtmaVm  to , size of [nQtma, nb]  
    end
    %QtmaBeqz num_G63
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dQtmaPertBeqz evaluated in x+pert
        [dSf_dQtma_PertBeqz, dSt_dQtma_PertBeqz] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertBeqz
        %2nd Derivatives of Sbr w.r.t. QtmaBeqz
        d2Sf_dBeqzQtma(:, k) = (dSf_dQtma_PertBeqz - dSf_dQtma).' * lam / pert;  %QtmaBeqz from, size of [nQtma, nBeqz]
        d2St_dBeqzQtma(:, k) = (dSt_dQtma_PertBeqz - dSt_dQtma).' * lam / pert;  %QtmaBeqz  to , size of [nQtma, nBeqz]
    end
    %QtmaBeqv num_G64
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dQtmaPertBeqv evaluated in x+pert
        [dSf_dQtma_PertBeqv, dSt_dQtma_PertBeqv] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertBeqv
        %2nd Derivatives of Sbr w.r.t. QtmaBeqv
        d2Sf_dBeqvQtma(:, k) = (dSf_dQtma_PertBeqv - dSf_dQtma).' * lam / pert;  %QtmaBeqv from, size of [nQtma, nBeqv]
        d2St_dBeqvQtma(:, k) = (dSt_dQtma_PertBeqv - dSt_dQtma).' * lam / pert;  %QtmaBeqv  to , size of [nQtma, nBeqv]
    end
    %QtmaPfsh num_G65
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dQtmaPertPfsh evaluated in x+pert
        [dSf_dQtma_PertPfsh, dSt_dQtma_PertPfsh] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertPfsh
        %2nd Derivatives of Sbr w.r.t. QtmaPfsh
        d2Sf_dPfshQtma(:, k) = (dSf_dQtma_PertPfsh - dSf_dQtma).' * lam / pert;  %QtmaPfsh from, size of [nQtma, nPfsh] 
        d2St_dPfshQtma(:, k) = (dSt_dQtma_PertPfsh - dSt_dQtma).' * lam / pert;  %QtmaPfsh  to , size of [nQtma, nPfsh] 
    end    
    
    %VxQtma num_G16 num_G26
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertQtma evaluated in x+pert
        [dSf_dV1_PertQtma, dSf_dV2_PertQtma, dSt_dV1_PertQtma, dSt_dV2_PertQtma, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        %2nd Derivatives of Sbr w.r.t. QtmaVx
        d2Sf_dQtmaVa(:, k) = (dSf_dV1_PertQtma - dSf_dV1).' * lam / pert;  %VaQtma from, size of [nb, nQtma] 
        d2Sf_dQtmaVm(:, k) = (dSf_dV2_PertQtma - dSf_dV2).' * lam / pert;  %VmQtma  to , size of [nb, nQtma] 
        d2St_dQtmaVa(:, k) = (dSt_dV1_PertQtma - dSt_dV1).' * lam / pert;  %VaQtma from, size of [nb, nQtma] 
        d2St_dQtmaVm(:, k) = (dSt_dV2_PertQtma - dSt_dV2).' * lam / pert;  %VmQtma  to , size of [nb, nQtma] 
    end
    %BeqzQtma num_G36
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertQtma evaluated in x+pert
        [dSf_dBeqz_PertQtma, dSt_dBeqz_PertQtma] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertQtma
        %2nd Derivatives of Sbr w.r.t. BeqzQtma
        d2Sf_dQtmaBeqz(:, k) = (dSf_dBeqz_PertQtma - dSf_dBeqz).' * lam / pert;  %BeqzQtma from, size of [nBeqz, nQtma] 
        d2St_dQtmaBeqz(:, k) = (dSt_dBeqz_PertQtma - dSt_dBeqz).' * lam / pert;  %BeqzQtma  to , size of [nBeqz, nQtma]  
    end
    %BeqvQtma num_G46
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqvPertQtma evaluated in x+pert
        [dSf_dBeqv_PertQtma, dSt_dBeqv_PertQtma] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertQtma
        %2nd Derivatives of Sbr w.r.t. BeqvQtma
        d2Sf_dQtmaBeqv(:, k) = (dSf_dBeqv_PertQtma - dSf_dBeqv).' * lam / pert;  %BeqvQtma from, size of [nBeqv, nQtma] 
        d2St_dQtmaBeqv(:, k) = (dSt_dBeqv_PertQtma - dSt_dBeqv).' * lam / pert;  %BeqvQtma  to , size of [nBeqv, nQtma]
    end
    %PfshQtma num_G56
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertQtma evaluated in x+pert
        [dSf_dPfsh_PertQtma, dSt_dPfsh_PertQtma] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertQtma
        %2nd Derivatives of Sbr w.r.t. PfshQtma   
        d2Sf_dQtmaPfsh(:, k) = (dSf_dPfsh_PertQtma - dSf_dPfsh).' * lam / pert;  %PfshQtma from, size of [nPfsh , nQtma] 
        d2St_dQtmaPfsh(:, k) = (dSt_dPfsh_PertQtma - dSt_dPfsh).' * lam / pert;  %PfshQtma  to , size of [nPfsh , nQtma] 
    end
    %Qtma2 num_G66
    for k=1:nQtma
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmaSel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dQtmaPertQtma evaluated in x+pert
        [dSf_dQtma_PertQtma, dSt_dQtma_PertQtma] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertQtma
        %2nd Derivatives of Sbr w.r.t. Qtma2   
        d2Sf_dQtma2(:, k) = (dSf_dQtma_PertQtma - dSf_dQtma).' * lam / pert;  %Qtma2 from, size of [nQtma , nQtma] 
        d2St_dQtma2(:, k) = (dSt_dQtma_PertQtma - dSt_dQtma).' * lam / pert;  %Qtma2  to , size of [nQtma , nQtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf16 = sparse(d2Sf_dQtmaVa);
num_Hf26 = sparse(d2Sf_dQtmaVm);
num_Hf36 = sparse(d2Sf_dQtmaBeqz);
num_Hf46 = sparse(d2Sf_dQtmaBeqv);
num_Hf56 = sparse(d2Sf_dQtmaPfsh);

num_Hf61 = sparse(d2Sf_dVaQtma);
num_Hf62 = sparse(d2Sf_dVmQtma);
num_Hf63 = sparse(d2Sf_dBeqzQtma);
num_Hf64 = sparse(d2Sf_dBeqvQtma);
num_Hf65 = sparse(d2Sf_dPfshQtma);

num_Hf66 = sparse(d2Sf_dQtma2);

num_Ht16 = sparse(d2St_dQtmaVa);
num_Ht26 = sparse(d2St_dQtmaVm);
num_Ht36 = sparse(d2St_dQtmaBeqz);
num_Ht46 = sparse(d2St_dQtmaBeqv);
num_Ht56 = sparse(d2St_dQtmaPfsh);

num_Ht61 = sparse(d2St_dVaQtma);
num_Ht62 = sparse(d2St_dVmQtma);
num_Ht63 = sparse(d2St_dBeqzQtma);
num_Ht64 = sparse(d2St_dBeqvQtma);
num_Ht65 = sparse(d2St_dPfshQtma);

num_Ht66 = sparse(d2St_dQtma2);



