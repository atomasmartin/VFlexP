function [num_dAf_dBeqx, num_dAt_dBeqx] = ...
                        dAbr_dBeqPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DABR_DBEQPERT  Partial derivatives of squared flow magnitudes w.r.t BEQ (Finite Differences Method).
%        [NUM_DF_DBEQX, NUM_DT_DBEQX] = ...
%                       DABR_DBEQPERT(BASEMVA, BUS, BRANCH, V, SIDE, PERT, VCART)
%   Returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   Beq using the Finite Differences Method.
%
%   Beq can be used either to control the Vdc to a certain set value Vfset 
%   or the Qf to match zero (zero constraint). So the derivatives are
%   separated for each function. The derivatives w.r.t. Beq will be 
%   chosen for either VSCI and VSCIIIz, or VSCII, or VSCI, VSCIIIz and VSCIII 
%   depending on the 3rd argument.
%   So that:
%
%   ctrl = 1 : Qf = 0,    Zero constraint   VSCI and VSCIIIz for Power Flow
%   ctrl = 2 : Vf = Vset, Vdc control       VSCII
%   ctrl = 3 : Qf = 0,    Zero constraint   VSCI ,VSCIIIz and VSCIII for OPF
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   The following explains the expressions used to form the matrices:
%
%   Let Af refer to the square of the apparent power at the "from" end of
%   each branch,
%
%       Af = abs(Sf).^2           %     At = abs(St).^2
%          = Sf .* conj(Sf)       %        = St .* conj(St)
%          = Pf.^2 + Qf.^2        %        = Pt.^2 + Qt.^2
%
%   then ...
%
%   Partial w.r.t Beqz 
%       dAf/dBeq = (AfPert - Af)/pert 
%       dAt/dBeq = (AtPert - At)/pert
%
%   Examples:
%       %% squared current magnitude
%       not coded yet
%
%       %% squared apparent power flow
%       [num_dAf_dBeqz, num_dAt_dBeqz] = ...
%                        dAbr_dBeqPert(baseMVA, bus, branch, V, 1, pert, vcart)
%       [num_dAf_dBeqv, num_dAt_dBeqv] = ...
%                        dAbr_dBeqPert(baseMVA, bus, branch, V, 2, pert, vcart)
%
%       %% squared real power flow
%       not coded yet
%
%   See also DIBR_DV, DSBR_DV, DABR_DV.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
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
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM

%% default input args
if nargin < 7
    vcart = 0;      %% default to polar coordinates
end
%% selection of VSC
if ctrl == 1 %VSC I and VSCIIIz
    iBeqx = find ( ( branch(:,CONV) == 1 | branch(:,CONV) == 3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC size[nBeqz,1]
elseif ctrl ==2 %VSC II
    iBeqx = find (branch(:,CONV) == ctrl & branch(:, BR_STATUS)==1 &  branch(:, VF_SET)~=0) ; %AAB- Find branch locations of VSC size[nBeqv,1]
elseif ctrl == 3 %VSC I, VSCIIIz and VSCIII
    iBeqx = find ( ( branch(:,CONV) == 1 | branch(:,CONV) == 3  | branch(:,CONV) == 4) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC size[nBeqz,1]
else
    error('dSbr_dBeq: VSC can only be control 1 (VSCI and VSCIIIz), 2 (VSCII), OR 3 (VSCI, VSCIIIz and VSCIII)')    
end  

%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nBeqx = length(iBeqx);      %% AAB- Number of VSC with active Beq

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dAbr_dBeqPert: Derivatives of Power balance equations w.r.t Beq using Finite Differences in cartasian has not been coded yet')    

else %AAB- Polar Version
   
    %Selector of active Beq 
    BeqAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    BeqAux(iBeqx) = 1; %AAB- Fill the selector with 1 where Beq is active
    diagBeqAux = sparse( diag(BeqAux) ); %AAB- Beq Selector [nl,nl]
    
    %Yf and Yt Original
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);     %AAB- obtain the Ybus, Yf, Yt
        
    %Sf and St evaluated in x
    Sf = diag(Cf*V) * conj(Yf * V);
    St = diag(Ct*V) * conj(Yt * V);
    
    %Af and At evaluated in x
    Af = Sf .* conj(Sf);
    At = St .* conj(St);    
     
    %Dimensionalize (Allocate for computational speed)
    num_dAf_dBeqx = sparse( zeros(nl,nBeqx) );
    num_dAt_dBeqx = sparse( zeros(nl,nBeqx) );
    
    for k=1:nBeqx
        PertSel=diagBeqAux(:,iBeqx(k)); %AAB- Selects the column of diagBeqsel representing only the active Beq
 
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.*PertSel); 
        
        %Yf and Yt Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert); %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sf and St evaluated in x+pert
        SfPert = diag(Cf*V) * conj(Yf_Pert * V);
        StPert = diag(Ct*V) * conj(Yt_Pert * V);
        
        %Af and At evaluated in x+pert
        AfPert = SfPert .* conj(SfPert);
        AtPert = StPert .* conj(StPert);  
         
        %Partials of Af and At w.r.t. Beq Finite differences f'(x) ~~ ( f(x+pert) - f(x) ) / pert 
        num_dAf_dBeqx(:, k) = (AfPert - Af )/ pert; %AAB- Final dAf_dBeq has a size of [nl, nBeqx]
        num_dAt_dBeqx(:, k) = (AtPert - At )/ pert; %AAB- Final dAt_dBeq has a size of [nl, nBeqx]
    end
end