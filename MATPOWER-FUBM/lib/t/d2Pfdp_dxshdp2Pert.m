function [num_Hf18, num_Hf28, num_Hf38, num_Hf48, num_Hf58, num_Hf68, num_Hf78, num_Hf81, num_Hf82, num_Hf83, num_Hf84, num_Hf85, num_Hf86, num_Hf87, num_Hf88] = d2Pfdp_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXSHDP2PERT  Computes 2nd derivatives of complex brch power flow "from" w.r.t. ShdpVa, ShdpVm, ShdpBeqz, ShdpBeqv, ShdpSh, Shdpqtma, Shdpvtma, VaShdp, VmShdp, BeqzShdp, BeqvShdp, ShShdp, qtmaShdp, vtmaShdp, ShdpShdp (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfShdpVa, HfShdpVm, HfShdpBz, HfShdpBv, HfShdpSh, HfShdpqtma, HfShdpvtma, HfVaShdp, HfVmShdp, HfBzShdp, HfBvShdp, HfShShdp, HfqtmaShdp, HfvtmaShdp, HfShdpShdp] = D2PFDP_DXSHDP2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 15 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector LAM with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HfShdpVa, HfShdpVm, HfShdpBz, HfShdpBv, HfShdpSh, HfShdpqtma, HfShdpvtma, HfVaShdp, HfVmShdp, HfBzShdp, HfBvShdp, HfShShdp, HfqtmaShdp, HfvtmaShdp, HfShdpShdp] = D2PFDP_DXSHDP2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfShdpVa, HfShdpVm, HfShdpBz, HfShdpBv, HfShdpSh, HfShdpqtma, HfShdpvtma, HfVaShdp, HfVmShdp, HfBzShdp, HfBvShdp, HfShShdp, HfqtmaShdp, HfvtmaShdp, HfShdpShdp] = D2PFDP_DXSHDP2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%                                                               
%       Here the output matrices correspond to:
%           HShdpVa   = d/dShdp (dPfdp_dVa.'   * mu)
%           HShdpVm   = d/dShdp (dPfdp_dVm.'   * mu)
%           HShdpBz   = d/dShdp (dPfdp_dBeqz.' * mu)
%           HShdpBv   = d/dShdp (dPfdp_dBeqv.' * mu)
%           HShdpSh   = d/dShdp (dPfdp_dSh.'   * mu)
%           HShdpqtma = d/dShdp (dPfdp_dqtma.' * mu)
%           HShdpvtma = d/dShdp (dPfdp_dvtma.' * mu)
%           HVaShdp   = d/dVa   (dPfdp_dShdp.' * mu)
%           HVmShdp   = d/dVm   (dPfdp_dShdp.' * mu)
%           HBzShdp   = d/dBeqz (dPfdp_dShdp.' * mu)
%           HBvShdp   = d/dBeqv (dPfdp_dShdp.' * mu)
%           HShShdp   = d/dSh   (dPfdp_dShdp.' * mu)
%           HqtmaShdp = d/dqtma (dPfdp_dShdp.' * mu)
%           HvtmaShdp = d/dvtma (dPfdp_dShdp.' * mu)
%           HShdpShdp = d/dShdp (dPfdp_dShdp.' * mu)
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
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 


%% Calculation of derivatives
if vcart
    error('d2Pfdp_dxshdp2Pert: Derivatives of Flow Limit equations w.r.t Theta_dp using Finite Differences Method in cartasian have not been coded yet')    

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
    [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4, vcart);    
    [dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3, vcart);
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    dPfdp_dBeqv = -real(dSf_dBeqv);
    dPfdp_dPfsh = -real(dSf_dPfsh);
    dPfdp_dQtma = -real(dSf_dQtma); 
    dPfdp_dVtma = -real(dSf_dVtma); 
    dPfdp_dPfdp = -real(dSf_dPfdp);     
    
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
    d2Pfdp_dPfdpVa   = sparse( zeros(nb,   nPfdp) ); 
    d2Pfdp_dPfdpVm   = sparse( zeros(nb,   nPfdp) );
    d2Pfdp_dPfdpBeqz = sparse( zeros(nBeqz,nPfdp) );
    d2Pfdp_dPfdpBeqv = sparse( zeros(nBeqv,nPfdp) );
    d2Pfdp_dPfdpPfsh = sparse( zeros(nPfsh,nPfdp) );   
    d2Pfdp_dPfdpQtma = sparse( zeros(nQtma,nPfdp) );  
    d2Pfdp_dPfdpVtma = sparse( zeros(nVtma,nPfdp) );
    
    d2Pfdp_dVaPfdp   = sparse( zeros(nPfdp,nb   ) ); 
    d2Pfdp_dVmPfdp   = sparse( zeros(nPfdp,nb   ) );
    d2Pfdp_dBeqzPfdp = sparse( zeros(nPfdp,nBeqz) );
    d2Pfdp_dBeqvPfdp = sparse( zeros(nPfdp,nBeqv) );
    d2Pfdp_dPfshPfdp = sparse( zeros(nPfdp,nPfsh) );   
    d2Pfdp_dQtmaPfdp = sparse( zeros(nPfdp,nQtma) );     
    d2Pfdp_dVtmaPfdp = sparse( zeros(nPfdp,nVtma) );
    
    d2Pfdp_dPfdp2    = sparse( zeros(nPfdp,nPfdp) );

    %PfdpVa num_G81
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dPfdp_PertVa, dSt_dPfdp_PertVa] = dSbr_dsh(branch, V1p, 3, vcart); %dSbr_dPfdpPertVa
        dPfdp_dPfdp_PertVa = -real(dSf_dPfdp_PertVa);         
        d2Pfdp_dVaPfdp(:, k) = (dPfdp_dPfdp_PertVa - dPfdp_dPfdp).' * lam / pert; %PfdpVa from, size of [nPfdp, nb] 
    end
    %PfdpVm num_G82
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dPfdp_PertVm, dSt_dPfdp_PertVm] = dSbr_dsh(branch, V2p, 3, vcart); %dSbr_dPfdpPertVm
        dPfdp_dPfdp_PertVm = -real(dSf_dPfdp_PertVm);        
        d2Pfdp_dVmPfdp(:, k) = (dPfdp_dPfdp_PertVm - dPfdp_dPfdp).' * lam / pert; %PfdpVm from, size of [nPfdp, nb] 
    end
    %PfdpBeqz num_G83
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVtmaPertBeqz evaluated in x+pert
        [dSf_dPfdp_PertBeqz, dSt_dPfdp_PertBeqz] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertBeqz
        dPfdp_dPfdp_PertBeqz = -real(dSf_dPfdp_PertBeqz);         
        %2nd Derivatives of Pfdp w.r.t. PfdpBeqz
        d2Pfdp_dBeqzPfdp(:, k) = (dPfdp_dPfdp_PertBeqz - dPfdp_dPfdp).' * lam / pert;  %PfdpBeqz from, size of [nPfdp, nBeqz]
    end
    %PfdpBeqv num_G84
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfdpPertBeqv evaluated in x+pert
        [dSf_dPfdp_PertBeqv, dSt_dPfdp_PertBeqv] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertBeqv
        dPfdp_dPfdp_PertBeqv = -real(dSf_dPfdp_PertBeqv);
        %2nd Derivatives of Pfdp w.r.t. PfdpBeqv
        d2Pfdp_dBeqvPfdp(:, k) = (dPfdp_dPfdp_PertBeqv - dPfdp_dPfdp).' * lam / pert;  %PfdpBeqv from, size of [nPfdp, nBeqv]
    end
    %PfdpPfsh num_G85
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfdpPertPfsh evaluated in x+pert
        [dSf_dPfdp_PertPfsh, dSt_dPfdp_PertPfsh] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertPfsh
        dPfdp_dPfdp_PertPfsh = -real(dSf_dPfdp_PertPfsh);        
        %2nd Derivatives of Pfdp w.r.t. PfdpPfsh
        d2Pfdp_dPfshPfdp(:, k) = (dPfdp_dPfdp_PertPfsh - dPfdp_dPfdp).' * lam / pert;  %PfdpPfsh from, size of [nPfdp, nPfsh] 
    end    
    %PfdpQtma num_G86
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfdpPertQtma evaluated in x+pert
        [dSf_dPfdp_PertQtma, dSt_dPfdp_PertQtma] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertQtma
        dPfdp_dPfdp_PertQtma = -real(dSf_dPfdp_PertQtma); 
        %2nd Derivatives of Pfdp w.r.t. PfdpQtma
        d2Pfdp_dQtmaPfdp(:, k) = (dPfdp_dPfdp_PertQtma - dPfdp_dPfdp).' * lam / pert;  %PfdpQtma from, size of [nPfdp, nQtma] 
    end 
    %PfdpVtma num_G87
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfdpPertVtma evaluated in x+pert
        [dSf_dPfdp_PertVtma, dSt_dPfdp_PertVtma] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertVtma
        dPfdp_dPfdp_PertVtma = -real(dSf_dPfdp_PertVtma); 
        %2nd Derivatives of Sbr w.r.t. PfdpVtma
        d2Pfdp_dVtmaPfdp(:, k) = (dPfdp_dPfdp_PertVtma - dPfdp_dPfdp).' * lam / pert;  %PfdpVtma from, size of [nPfdp, nVtma] 
    end 
    
    %VxVtma num_G18 num_G28
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertPfdp evaluated in x+pert
        [dSf_dV1_PertPfdp, dSf_dV2_PertPfdp, dSt_dV1_PertPfdp dSt_dV2_PertPfdp, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        
        %dPfdp_dVxPertVtma evaluated in x+pert
        dPfdp_dV1_PertPfdp = -real(dSf_dV1_PertPfdp);
        dPfdp_dV2_PertPfdp = -real(dSf_dV2_PertPfdp) + Kdp.*( dVmf_dVm );
                
        %2nd Derivatives of Pfdp w.r.t. VtmaVx
        d2Pfdp_dPfdpVa(:, k) = (dPfdp_dV1_PertPfdp - dPfdp_dV1).' * lam / pert;  %VaVtma from, size of [nb, nPfdp] 
        d2Pfdp_dPfdpVm(:, k) = (dPfdp_dV2_PertPfdp - dPfdp_dV2).' * lam / pert;  %VmVtma  to , size of [nb, nPfdp] 
    end
    %BeqzPfdp num_G38
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertPfdp evaluated in x+pert
        [dSf_dBeqz_PertPfdp, dSt_dBeqz_PertPfdp] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertPfdp
        dPfdp_dBeqz_PertPfdp = -real(dSf_dBeqz_PertPfdp);
        %2nd Derivatives of Sbr w.r.t. BeqzPfdp
        d2Pfdp_dPfdpBeqz(:, k) = (dPfdp_dBeqz_PertPfdp - dPfdp_dBeqz).' * lam / pert;  %BeqzPfdp from, size of [nBeqz, nPfdp] 
    end
    %BeqvPfdp num_G48
    for k=1:nPfdp 
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqvPertPfdp evaluated in x+pert
        [dSf_dBeqv_PertPfdp, dSt_dBeqv_PertPfdp] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertPfdp
        dPfdp_dBeqv_PertPfdp = -real(dSf_dBeqv_PertPfdp);
        %2nd Derivatives of dPfdp_dBeqv_PertVtma = -real(dSf_dBeqv_PertVtma); w.r.t. BeqvPfdp
        d2Pfdp_dPfdpBeqv(:, k) = (dPfdp_dBeqv_PertPfdp - dPfdp_dBeqv).' * lam / pert;  %BeqvPfdp from, size of [nBeqv, nPfdp] 
    end
    %PfshPfdp num_G57
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertPfdp evaluated in x+pert
        [dSf_dPfsh_PertPfdp, dSt_dPfsh_PertPfdp] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertPfdp
        dPfdp_dPfsh_PertPfdp = -real(dSf_dPfsh_PertPfdp);
        %2nd Derivatives of Pfdp w.r.t. PfshPfdp   
        d2Pfdp_dPfdpPfsh(:, k) = (dPfdp_dPfsh_PertPfdp - dPfdp_dPfsh).' * lam / pert;  %PfshPfdp from, size of [nPfsh , nPfdp] 
    end
    %QtmaPfdp num_G68
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %dSbr_dQtmaPertPfdp evaluated in x+pert
        [dSf_dQtma_PertPfdp, dSt_dQtma_PertPfdp] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertPfdp
        dPfdp_dQtma_PertPfdp = -real(dSf_dQtma_PertPfdp);
        %2nd Derivatives of Sbr w.r.t. QtmaPfdp   
        d2Pfdp_dPfdpQtma(:, k) = (dPfdp_dQtma_PertPfdp - dPfdp_dQtma).' * lam / pert;  %QtmaPfdp from, size of [nQtma , nPfdp] 
    end
    %VtmaPfdp num_G78
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %dSbr_dVtmaPertPfdp evaluated in x+pert
        [dSf_dVtma_PertPfdp, dSt_dVtma_PertPfdp] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertPfdp
        dPfdp_dVtma_PertPfdp = -real(dSf_dVtma_PertPfdp);
        %2nd Derivatives of Pfdp w.r.t. VtmaPfdp   
        d2Pfdp_dPfdpVtma(:, k) = (dPfdp_dVtma_PertPfdp - dPfdp_dVtma).' * lam / pert;  %VtmaPfdp from, size of [nVtma , nPfdp] 
    end
    %Pfdp2 num_G88
    for k=1:nPfdp
        PertSel=diagShdpSel(:,iPfdp(k)); %AAB- Selects the column of diagShdpsel representing only the active Shdp
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_dp in the Perturbed branch (One Pfdp at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);        
        %dSbr_dPfdpPertPfdp evaluated in x+pert
        [dSf_dPfdp_PertPfdp, dSt_dPfdp_PertPfdp] = dSbr_dsh(branch_Pert, V, 3, vcart); %dSbr_dPfdpPertPfdp
        dPfdp_dPfdp_PertPfdp = -real(dSf_dPfdp_PertPfdp);
        %2nd Derivatives of Sbr w.r.t. Pfdp2   
        d2Pfdp_dPfdp2(:, k) = (dPfdp_dPfdp_PertPfdp - dPfdp_dPfdp).' * lam / pert;  %Pfdp2 from, size of [nPfdp , nPfdp] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf18 = sparse(d2Pfdp_dPfdpVa);
num_Hf28 = sparse(d2Pfdp_dPfdpVm);
num_Hf38 = sparse(d2Pfdp_dPfdpBeqz);
num_Hf48 = sparse(d2Pfdp_dPfdpBeqv);
num_Hf58 = sparse(d2Pfdp_dPfdpPfsh);
num_Hf68 = sparse(d2Pfdp_dPfdpQtma);
num_Hf78 = sparse(d2Pfdp_dPfdpVtma);

num_Hf81 = sparse(d2Pfdp_dVaPfdp);
num_Hf82 = sparse(d2Pfdp_dVmPfdp);
num_Hf83 = sparse(d2Pfdp_dBeqzPfdp);
num_Hf84 = sparse(d2Pfdp_dBeqvPfdp);
num_Hf85 = sparse(d2Pfdp_dPfshPfdp);
num_Hf86 = sparse(d2Pfdp_dQtmaPfdp);
num_Hf87 = sparse(d2Pfdp_dVtmaPfdp);

num_Hf88 = sparse(d2Pfdp_dPfdp2);
