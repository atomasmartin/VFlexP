function [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55] = d2Pfdp_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXSH2PERT  Computes 2nd derivatives of Droop Control "from" w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh] = D2PFDP_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 9 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq and Sh of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh] = D2PFDP_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh] = d2PFDP_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%
%       Here the output matrices correspond to:
%           HShVa = d/dSh   (dPfdp_dVa.'   * mu)
%           HShVm = d/dSh   (dPfdp_dVm.'   * mu)
%           HShBz = d/dSh   (dPfdp_dBeqz.' * mu)
%           HShBv = d/dSh   (dPfdp_dBeqv.' * mu)
%           HVaSh = d/dVa   (dPfdp_dSh.'   * mu)
%           HVmSh = d/dVm   (dPfdp_dSh.'   * mu)
%           HBzSh = d/dBeqz (dPfdp_dSh.'   * mu)
%           HBvSh = d/dBeqv (dPfdp_dSh.'   * mu)
%           HShSh = d/dSh   (dPfdp_dSh.'   * mu)
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
    ALPH1, ALPH2, ALPH3 KDP] = idx_brch;%<<AAB-extra fields for FUBM
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

%% identifier of AC/DC grids
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)& (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

%% Calculation of derivatives
if vcart
    error('d2Pfdp_dxsh2Pert: Derivatives of Flow Limit equations w.r.t Theta_sh using Finite Differences Method in cartasian have not been coded yet')    

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
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    dPfdp_dBeqv = -real(dSf_dBeqv);
    dPfdp_dPfsh = -real(dSf_dPfsh);
    
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
      
    %Selector of active Beqv 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqvAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
        
    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
   
    %Dimensionalize (Allocate for computational speed)  
    d2Pfdp_dPfshVa   = sparse( zeros(nb,   nPfsh) ); 
    d2Pfdp_dPfshVm   = sparse( zeros(nb,   nPfsh) );
    d2Pfdp_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2Pfdp_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );

    d2Pfdp_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2Pfdp_dVmPfsh   = sparse( zeros(nPfsh,nb   ) );
    d2Pfdp_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2Pfdp_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
    
    d2Pfdp_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );

    %PfShVa num_G51
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dPfsh_PertVa, dSt_dPfsh_PertVa] = dSbr_dsh(branch, V1p, 1, vcart); %dSbr_dPfshPertVa
        dPfdp_dPfsh_PertVa = -real(dSf_dPfsh_PertVa);     
        d2Pfdp_dVaPfsh(:, k) = (dPfdp_dPfsh_PertVa - dPfdp_dPfsh).' * lam / pert; %shVa from, size of [nPfsh, nb]
    end
    %PfshVm num_G52
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dPfsh_PertVm, dSt_dPfsh_PertVm] = dSbr_dsh(branch, V2p, 1, vcart); %dSbr_dPfshPertVm
        dPfdp_dPfsh_PertVm = -real(dSf_dPfsh_PertVm);     
        d2Pfdp_dVmPfsh(:, k) = (dPfdp_dPfsh_PertVm - dPfdp_dPfsh).' * lam / pert; %shVm from, size of [nPfsh, nb]
    end
    %PfshBeqz num_G53
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertBeqz evaluated in x+pert
        [dSf_dPfsh_PertBeqz, dSt_dPfsh_PertBeqz] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertBeqz
        dPfdp_dPfsh_PertBeqz = -real(dSf_dPfsh_PertBeqz);         
        %2nd Derivatives of Pfdp w.r.t. PfshBeqz
        d2Pfdp_dBeqzPfsh(:, k) = (dPfdp_dPfsh_PertBeqz - dPfdp_dPfsh).' * lam / pert;  %PfshBeqz from, size of [nPfsh, nBeqz]
    end
    %PfshBeqv num_G54
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertBeqv evaluated in x+pert
        [dSf_dPfsh_PertBeqv,dSt_dPfsh_PertBeqv] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertBeqv
        dPfdp_dPfsh_PertBeqv = -real(dSf_dPfsh_PertBeqv);         
        %2nd Derivatives of Sbr w.r.t. PfshBeqv
        d2Pfdp_dBeqvPfsh(:, k) = (dPfdp_dPfsh_PertBeqv - dPfdp_dPfsh).' * lam / pert;  %PfshBeqv from, size of [nPfsh, nBeqv]
    end
    %VxPfsh num_G15 num_G25
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertPfsh evaluated in x+pert
        [dSf_dV1_PertPfsh, dSf_dV2_PertPfsh, dSt_dV1_PertPfsh, dSt_dV2_PertPfsh, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        
        %dPfdp_dVxPertPfsh evaluated in x+pert
        dPfdp_dV1_PertPfsh = -real(dSf_dV1_PertPfsh);
        dPfdp_dV2_PertPfsh = -real(dSf_dV2_PertPfsh) + Kdp.*( dVmf_dVm );
        
        %2nd Derivatives of Pfdp w.r.t. PfshVx
        d2Pfdp_dPfshVa(:, k) = (dPfdp_dV1_PertPfsh - dPfdp_dV1).' * lam / pert;  %VaPfsh from, size of [nb, nPfsh] 
        d2Pfdp_dPfshVm(:, k) = (dPfdp_dV2_PertPfsh - dPfdp_dV2).' * lam / pert;  %VmPfsh  to , size of [nb, nPfsh] 
    end
    %BeqzPfsh num_G35
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertPfsh evaluated in x+pert
        [dSf_dBeqz_PertPfsh, dSt_dBeqz_PertPfsh] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertPfsh
        dPfdp_dBeqz_PertPfsh = -real(dSf_dBeqz_PertPfsh);       
        %2nd Derivatives of Sbr w.r.t. BeqzPfsh
        d2Pfdp_dPfshBeqz(:, k) = (dPfdp_dBeqz_PertPfsh - dPfdp_dBeqz).' * lam / pert;  %BeqzPfsh from, size of [nBeqz, nPfsh] 
    end
    %BeqvPfsh num_G45
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqvPertPfsh evaluated in x+pert
        [dSf_dBeqv_PertPfsh, dSt_dBeqv_PertPfsh] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertPfsh
        dPfdp_dBeqv_PertPfsh = -real(dSf_dBeqv_PertPfsh);  
        
        %2nd Derivatives of Pfdp w.r.t. BeqvPfsh
        d2Pfdp_dPfshBeqv(:, k) = (dPfdp_dBeqv_PertPfsh - dPfdp_dBeqv).' * lam / pert;  %BeqvPfsh from, size of [nBeqv, nPfsh] 
    end
    %PfshPfsh num_G55
    for k=1:nPfsh
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShSel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertPfsh evaluated in x+pert
        [dSf_dPfsh_PertPfsh, dSt_dPfsh_PertPfsh] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertPfsh
        dPfdp_dPfsh_PertPfsh = -real(dSf_dPfsh_PertPfsh);        
        
        %2nd Derivatives of Pfdp w.r.t. Pfsh2   
        d2Pfdp_dPfsh2(:, k) = (dPfdp_dPfsh_PertPfsh - dPfdp_dPfsh).' * lam / pert;  %PfshPfsh from, size of [nPfsh , nPfsh] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf15 = sparse(d2Pfdp_dPfshVa);
num_Hf25 = sparse(d2Pfdp_dPfshVm);
num_Hf35 = sparse(d2Pfdp_dPfshBeqz);
num_Hf45 = sparse(d2Pfdp_dPfshBeqv);

num_Hf51 = sparse(d2Pfdp_dVaPfsh);
num_Hf52 = sparse(d2Pfdp_dVmPfsh);
num_Hf53 = sparse(d2Pfdp_dBeqzPfsh);
num_Hf54 = sparse(d2Pfdp_dBeqvPfsh);

num_Hf55 = sparse(d2Pfdp_dPfsh2);
