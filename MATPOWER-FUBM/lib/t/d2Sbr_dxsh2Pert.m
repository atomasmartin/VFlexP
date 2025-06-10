function [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55,...
          num_Ht15, num_Ht25, num_Ht35, num_Ht45, num_Ht51, num_Ht52, num_Ht53, num_Ht54, num_Ht55] = d2Sbr_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2SF_DXSH2PERT  Computes 2nd derivatives of complex brch power flow "from" w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh,...
%    HtShVa, HtShVm, HtShBz, HtShBv, HtVaSh, HtVmSh, HtBzSh, HtBvSh, HtShSh] = D2SF_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 18 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq and Sh of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh,...
%    HtShVa, HtShVm, HtShBz, HtShBv, HtVaSh, HtVmSh, HtBzSh, HtBvSh, HtShSh] = D2SF_DXSH2PERT(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfShVa, HfShVm, HfShBz, HfShBv, HfVaSh, HfVmSh, HfBzSh, HfBvSh, HfShSh,...
%    HtShVa, HtShVm, HtShBz, HtShBv, HtVaSh, HtVmSh, HtBzSh, HtBvSh, HtShSh] = d2Sbr_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%
%       Here the output matrices correspond to:
%           HShVa = d/dSh   (dSf_dVa.'   * mu)
%           HShVm = d/dSh   (dSf_dVm.'   * mu)
%           HShBz = d/dSh   (dSf_dBeqz.' * mu)
%           HShBv = d/dSh   (dSf_dBeqv.' * mu)
%           HVaSh = d/dVa   (dSf_dSh.'   * mu)
%           HVmSh = d/dVm   (dSf_dSh.'   * mu)
%           HBzSh = d/dBeqz (dSf_dSh.'   * mu)
%           HBvSh = d/dBeqv (dSf_dSh.'   * mu)
%           HShSh = d/dSh   (dSf_dSh.'   * mu)
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

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift

%% Calculation of derivatives
if vcart
    error('d2Sf_dxsh2Pert: Derivatives of Flow Limit equations w.r.t Theta_sh using Finite Differences Method in cartasian have not been coded yet')    

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
   
    %Dimensionalize (Allocate for computational speed)  
    d2Sf_dPfshVa   = sparse( zeros(nb,   nPfsh) ); 
    d2Sf_dPfshVm   = sparse( zeros(nb,   nPfsh) );
    d2Sf_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2Sf_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );

    d2Sf_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2Sf_dVmPfsh   = sparse( zeros(nPfsh,nb   ) );
    d2Sf_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2Sf_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
    
    d2Sf_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );

    d2St_dPfshVa   = sparse( zeros(nb,   nPfsh) ); 
    d2St_dPfshVm   = sparse( zeros(nb,   nPfsh) );
    d2St_dPfshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2St_dPfshBeqv = sparse( zeros(nBeqv,nPfsh) );

    d2St_dVaPfsh   = sparse( zeros(nPfsh,nb   ) ); 
    d2St_dVmPfsh   = sparse( zeros(nPfsh,nb   ) );
    d2St_dBeqzPfsh = sparse( zeros(nPfsh,nBeqz) );
    d2St_dBeqvPfsh = sparse( zeros(nPfsh,nBeqv) );
    
    d2St_dPfsh2    = sparse( zeros(nPfsh,nPfsh) );    
    
    %PfShVa num_G51
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dPfsh_PertVa, dSt_dPfsh_PertVa] = dSbr_dsh(branch, V1p, 1, vcart); %dSbr_dPfshPertVa
        d2Sf_dVaPfsh(:, k) = (dSf_dPfsh_PertVa - dSf_dPfsh).' * lam / pert; %shVa from, size of [nPfsh, nb]
        d2St_dVaPfsh(:, k) = (dSt_dPfsh_PertVa - dSt_dPfsh).' * lam / pert; %shVa  to , size of [nPfsh, nb]
    end
    %PfshVm num_G52
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dPfsh_PertVm, dSt_dPfsh_PertVm] = dSbr_dsh(branch, V2p, 1, vcart); %dSbr_dPfshPertVm
        d2Sf_dVmPfsh(:, k) = (dSf_dPfsh_PertVm - dSf_dPfsh).' * lam / pert; %shVm from, size of [nPfsh, nb]
        d2St_dVmPfsh(:, k) = (dSt_dPfsh_PertVm - dSt_dPfsh).' * lam / pert; %shVm  to , size of [nPfsh, nb]
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
        %2nd Derivatives of Sbr w.r.t. PfshBeqz
        d2Sf_dBeqzPfsh(:, k) = (dSf_dPfsh_PertBeqz - dSf_dPfsh).' * lam / pert;  %PfshBeqz from, size of [nPfsh, nBeqz]
        d2St_dBeqzPfsh(:, k) = (dSt_dPfsh_PertBeqz - dSt_dPfsh).' * lam / pert;  %PfshBeqz  to , size of [nPfsh, nBeqz]
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
        %2nd Derivatives of Sbr w.r.t. PfshBeqv
        d2Sf_dBeqvPfsh(:, k) = (dSf_dPfsh_PertBeqv - dSf_dPfsh).' * lam / pert;  %PfshBeqv from, size of [nPfsh, nBeqv]
        d2St_dBeqvPfsh(:, k) = (dSt_dPfsh_PertBeqv - dSt_dPfsh).' * lam / pert;  %PfshBeqv  to , size of [nPfsh, nBeqv]
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
        %2nd Derivatives of Sbr w.r.t. PfshVx
        d2Sf_dPfshVa(:, k) = (dSf_dV1_PertPfsh - dSf_dV1).' * lam / pert;  %VaPfsh from, size of [nb, nPfsh] 
        d2Sf_dPfshVm(:, k) = (dSf_dV2_PertPfsh - dSf_dV2).' * lam / pert;  %VmPfsh  to , size of [nb, nPfsh] 
        d2St_dPfshVa(:, k) = (dSt_dV1_PertPfsh - dSt_dV1).' * lam / pert;  %VaPfsh from, size of [nb, nPfsh] 
        d2St_dPfshVm(:, k) = (dSt_dV2_PertPfsh - dSt_dV2).' * lam / pert;  %VmPfsh  to , size of [nb, nPfsh]         
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
        %2nd Derivatives of Sbr w.r.t. BeqzPfsh
        d2Sf_dPfshBeqz(:, k) = (dSf_dBeqz_PertPfsh - dSf_dBeqz).' * lam / pert;  %BeqzPfsh from, size of [nBeqz, nPfsh] 
        d2St_dPfshBeqz(:, k) = (dSt_dBeqz_PertPfsh - dSt_dBeqz).' * lam / pert;  %BeqzPfsh  to , size of [nBeqz, nPfsh] 
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
        %2nd Derivatives of Sbr w.r.t. BeqvPfsh
        d2Sf_dPfshBeqv(:, k) = (dSf_dBeqv_PertPfsh - dSf_dBeqv).' * lam / pert;  %BeqvPfsh from, size of [nBeqv, nPfsh] 
        d2St_dPfshBeqv(:, k) = (dSt_dBeqv_PertPfsh - dSt_dBeqv).' * lam / pert;  %BeqvPfsh  to , size of [nBeqv, nPfsh] 
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
        %2nd Derivatives of Sbus w.r.t. Pfsh2   
        d2Sf_dPfsh2(:, k) = (dSf_dPfsh_PertPfsh - dSf_dPfsh).' * lam / pert;  %PfshPfsh from, size of [nPfsh , nPfsh] 
        d2St_dPfsh2(:, k) = (dSt_dPfsh_PertPfsh - dSt_dPfsh).' * lam / pert;  %PfshPfsh  to , size of [nPfsh , nPfsh]
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf15 = sparse(d2Sf_dPfshVa);
num_Hf25 = sparse(d2Sf_dPfshVm);
num_Hf35 = sparse(d2Sf_dPfshBeqz);
num_Hf45 = sparse(d2Sf_dPfshBeqv);

num_Hf51 = sparse(d2Sf_dVaPfsh);
num_Hf52 = sparse(d2Sf_dVmPfsh);
num_Hf53 = sparse(d2Sf_dBeqzPfsh);
num_Hf54 = sparse(d2Sf_dBeqvPfsh);

num_Hf55 = sparse(d2Sf_dPfsh2);

num_Ht15 = sparse(d2St_dPfshVa);
num_Ht25 = sparse(d2St_dPfshVm);
num_Ht35 = sparse(d2St_dPfshBeqz);
num_Ht45 = sparse(d2St_dPfshBeqv);

num_Ht51 = sparse(d2St_dVaPfsh);
num_Ht52 = sparse(d2St_dVmPfsh);
num_Ht53 = sparse(d2St_dBeqzPfsh);
num_Ht54 = sparse(d2St_dBeqvPfsh);

num_Ht55 = sparse(d2St_dPfsh2);









