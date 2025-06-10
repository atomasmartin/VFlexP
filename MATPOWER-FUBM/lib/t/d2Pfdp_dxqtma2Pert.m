function [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66] = d2Pfdp_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXQTMA2PERT  Computes 2nd derivatives of Droop Control "from" w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma] = D2PFDP_DXQTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 11 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector MU with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma] = D2PFDP_DXQTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfqtmaVa, HfqtmaVm, HfqtmaBz, HfqtmaBv, HfqtmaSh, HfVaqtma, HfVmqtma, HfBzqtma, HfBvqtma, HfShqtma, Hfqtmaqtma] = d2Pfdp_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%
%       Here the output matrices correspond to:
%           HqtmaVa   = d/dqtma (dPfdp_dVa.'   * mu)
%           HqtmaVm   = d/dqtma (dPfdp_dVm.'   * mu)
%           HqtmaBz   = d/dqtma (dPfdp_dBeqz.' * mu)
%           HqtmaBv   = d/dqtma (dPfdp_dBeqv.' * mu)
%           HqtmaSh   = d/dqtma (dPfdp_dSh.'   * mu)
%           HVaqtma   = d/dVa   (dPfdp_dqtma.' * mu)
%           HVmqtma   = d/dVm   (dPfdp_dqtma.' * mu)
%           HBzqtma   = d/dBeqz (dPfdp_dqtma.' * mu)
%           HBvqtma   = d/dBeqv (dPfdp_dqtma.' * mu)
%           HShqtma   = d/dSh   (dPfdp_dqtma.' * mu)
%           Hqtmaqtma = d/dqtma (dPfdp_dqtma.' * mu)
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
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)& (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

%% Calculation of derivatives
if vcart
    error('d2Pfdp_dxma2Pert: Derivatives of Flow Limit equations w.r.t ma using Finite Differences Method in cartasian have not been coded yet')    

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
    [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1, vcart);
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, vcart);
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    dPfdp_dBeqv = -real(dSf_dBeqv);
    dPfdp_dPfsh = -real(dSf_dPfsh);
    dPfdp_dQtma = -real(dSf_dQtma);    
   
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
    d2Pfdp_dQtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Pfdp_dQtmaVm   = sparse( zeros(nb,   nQtma) );
    d2Pfdp_dQtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Pfdp_dQtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Pfdp_dQtmaPfsh = sparse( zeros(nPfsh,nQtma) ); 

    d2Pfdp_dVaQtma   = sparse( zeros(nQtma,nb   ) ); 
    d2Pfdp_dVmQtma   = sparse( zeros(nQtma,nb   ) );
    d2Pfdp_dBeqzQtma = sparse( zeros(nQtma,nBeqz) );
    d2Pfdp_dBeqvQtma = sparse( zeros(nQtma,nBeqv) );
    d2Pfdp_dPfshQtma = sparse( zeros(nQtma,nPfsh) ); 
    
    d2Pfdp_dQtma2    = sparse( zeros(nQtma,nQtma) );

    %QtmaVa num_G61
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dQtma_PertVa, dSt_dQtma_PertVa] = dSbr_dma(branch, V1p, 2, vcart); %dSbr_dQtmaPertVa
        dPfdp_dQtma_PertVa = -real(dSf_dQtma_PertVa);     
        d2Pfdp_dVaQtma(:, k) = (dPfdp_dQtma_PertVa - dPfdp_dQtma).' * lam / pert; %QtmaVa from, size of [nQtma, nb] 
    end
    %QtmaVm num_G62
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dQtma_PertVm, dSt_dQtma_PertVm] = dSbr_dma(branch, V2p, 2, vcart); %dSbr_dQtmaPertVm
        dPfdp_dQtma_PertVm = -real(dSf_dQtma_PertVm);     
        d2Pfdp_dVmQtma(:, k) = (dPfdp_dQtma_PertVm - dPfdp_dQtma).' * lam / pert; %QtmaVm from, size of [nQtma, nb] 
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
        dPfdp_dQtma_PertBeqz = -real(dSf_dQtma_PertBeqz);                 
        %2nd Derivatives of Pfdp w.r.t. QtmaBeqz
        d2Pfdp_dBeqzQtma(:, k) = (dPfdp_dQtma_PertBeqz - dPfdp_dQtma).' * lam / pert;  %QtmaBeqz from, size of [nQtma, nBeqz]
        d2Pfdp_dBeqzQtma(:, k) = (dPfdp_dQtma_PertBeqz - dPfdp_dQtma).' * lam / pert;  %QtmaBeqz  to , size of [nQtma, nBeqz]
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
        dPfdp_dQtma_PertBeqv = -real(dSf_dQtma_PertBeqv);         
        %2nd Derivatives of Sbr w.r.t. QtmaBeqv
        d2Pfdp_dBeqvQtma(:, k) = (dPfdp_dQtma_PertBeqv - dPfdp_dQtma).' * lam / pert;  %QtmaBeqv from, size of [nQtma, nBeqv]
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
        dPfdp_dQtma_PertPfsh = -real(dSf_dQtma_PertPfsh);               
        %2nd Derivatives of Sbr w.r.t. QtmaPfsh
        d2Pfdp_dPfshQtma(:, k) = (dPfdp_dQtma_PertPfsh - dPfdp_dQtma).' * lam / pert;  %QtmaPfsh from, size of [nQtma, nPfsh] 
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

        %dPfdp_dVxPertQtma evaluated in x+pert
        dPfdp_dV1_PertQtma = -real(dSf_dV1_PertQtma);
        dPfdp_dV2_PertQtma = -real(dSf_dV2_PertQtma) + Kdp.*( dVmf_dVm );
        
        %2nd Derivatives of Pfdp w.r.t. QtmaVx
        d2Pfdp_dQtmaVa(:, k) = (dPfdp_dV1_PertQtma - dPfdp_dV1).' * lam / pert;  %VaQtma from, size of [nb, nQtma] 
        d2Pfdp_dQtmaVm(:, k) = (dPfdp_dV2_PertQtma - dPfdp_dV2).' * lam / pert;  %VmQtma  to , size of [nb, nQtma] 
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
        dPfdp_dBeqz_PertQtma = -real(dSf_dBeqz_PertQtma); 
        %2nd Derivatives of Sbr w.r.t. BeqzQtma
        d2Pfdp_dQtmaBeqz(:, k) = (dPfdp_dBeqz_PertQtma - dPfdp_dBeqz).' * lam / pert;  %BeqzQtma from, size of [nBeqz, nQtma] 
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
        dPfdp_dBeqv_PertQtma = -real(dSf_dBeqv_PertQtma);
        %2nd Derivatives of Sbr w.r.t. BeqvQtma
        d2Pfdp_dQtmaBeqv(:, k) = (dPfdp_dBeqv_PertQtma - dPfdp_dBeqv).' * lam / pert;  %BeqvQtma from, size of [nBeqv, nQtma] 
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
        dPfdp_dPfsh_PertQtma = -real(dSf_dPfsh_PertQtma);
        %2nd Derivatives of Pfdp w.r.t. PfshQtma   
        d2Pfdp_dQtmaPfsh(:, k) = (dPfdp_dPfsh_PertQtma - dPfdp_dPfsh).' * lam / pert;  %PfshQtma from, size of [nPfsh , nQtma] 
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
        dPfdp_dQtma_PertQtma = -real(dSf_dQtma_PertQtma);
        %2nd Derivatives of Sbr w.r.t. Qtma2   
        d2Pfdp_dQtma2(:, k) = (dPfdp_dQtma_PertQtma - dPfdp_dQtma).' * lam / pert;  %Qtma2 from, size of [nQtma , nQtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf16 = sparse(d2Pfdp_dQtmaVa);
num_Hf26 = sparse(d2Pfdp_dQtmaVm);
num_Hf36 = sparse(d2Pfdp_dQtmaBeqz);
num_Hf46 = sparse(d2Pfdp_dQtmaBeqv);
num_Hf56 = sparse(d2Pfdp_dQtmaPfsh);

num_Hf61 = sparse(d2Pfdp_dVaQtma);
num_Hf62 = sparse(d2Pfdp_dVmQtma);
num_Hf63 = sparse(d2Pfdp_dBeqzQtma);
num_Hf64 = sparse(d2Pfdp_dBeqvQtma);
num_Hf65 = sparse(d2Pfdp_dPfshQtma);

num_Hf66 = sparse(d2Pfdp_dQtma2);
