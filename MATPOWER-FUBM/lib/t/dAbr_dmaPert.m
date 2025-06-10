function [num_dAf_dXxma, num_dAt_dXxma] = ...
                        dAbr_dmaPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DABR_DMAPERT  Partial derivatives of squared flow magnitudes w.r.t ma/tap (Finite Differences Method).
%        [NUM_DAF_DXXMA, NUM_DAT_DXXMA] = ...
%                       DABR_DMAPERT(DFF_DXXMA, DFT_DXXMA, FF, FT)
%   Returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   ma/tap using the Finite Differences Method
%
%   ma can be used either to control the Qf, Qt, Vf or Vt to a certain set 
%   value Qfset, Qtset, Vfset or Vtset respectively. Thus, the derivatives 
%   are separated for each function. The derivatives w.r.t. ma for a
%   certain control will be chosen depending on the 5rd argument "ctrl".
%   So that:
%
%   ctrl = 1 : Qf = Qfset, "from side", Transformers (not for VSC)
%   ctrl = 2 : Qt = Qtset, "to side",   Transformers and VSC
%   ctrl = 3 : Vf = Vfset, "from side", Transformers (not for VSC)
%   ctrl = 4 : Vt = Vtset, "to side",   Transformers and VSC
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
%   Partial w.r.t Theta_sh 
%       dAf/dma = (AfPert - Af)/pert 
%       dAt/dma = (AtPert - At)/pert
%
%   Examples:
%       %% squared current magnitude
%       not coded yet
%
%       %% squared apparent power flow
%       [num_dAf_dQtma, num_dAt_dQtma] = ...
%                        dAbr_dmaPert(baseMVA, bus, branch, V, 2, pert, vcart)
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
    ALPH1, ALPH2, ALPH3] = idx_brch;%<<AAB-extra fields for FUBM

%% default input args
if nargin < 7
    vcart = 0;      %% default to polar coordinates
end
%% control selection
if ctrl == 1     %Qf
    iXxma = find (branch(:,QF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VF_SET)==0 & branch(:, CONV)==0  ); %AAB- Find branch locations of Qf control size[nQfma,1] %Transformers
elseif ctrl == 2 %Qt
    iXxma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations of Qt control size[nQtsh,1] %Transformers and VSC
elseif ctrl == 3 %Vf
    iXxma = find (branch(:,VF_SET)~=0 & branch(:, BR_STATUS)==1 & branch(:, CONV)==0 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) ); %AAB- Find branch locations of Qt control size[nQtsh,1] %Transformers
elseif ctrl == 4 %Vt
    iXxma = find (branch(:,VT_SET)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) ); %AAB- Find branch locations of Qt control size[nQtsh,1] %Transformers and VSC
else
    error('dSbus_dsh: Control type can only be type 1 (Qf), 2 (Qt), 3(Vf), or 4(Vt)')    
end  
%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nXxma = length(iXxma);      %% AAB- Number of elements with Voltage or Reactive power controlled by ma  

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dAbr_dma: Derivatives of Power balance equations w.r.t ma using Finite Differences in cartasian have not been coded yet')    

else %AAB- Polar Version
    
    %Selector of active ma for the specified control 
    maAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    maAux(iXxma) = 1; %AAB- Fill the selector with "1" where ma is active
    diagmaAux = sparse( diag(maAux) );  %AAB- Diagonal of ma selector  size [nl,nl]
    
    %Yf and Yt Original
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);     %AAB- obtain the Ybus, Yf, Yt
        
    %Sf and St evaluated in x
    Sf = diag(Cf*V) * conj(Yf * V);
    St = diag(Ct*V) * conj(Yt * V);
    
    %Af and At evaluated in x
    Af = Sf .* conj(Sf);
    At = St .* conj(St);    
            
    %Dimensionalize (Allocate for computational speed)
    num_dAf_dXxma = sparse( zeros(nl,nXxma) );
    num_dAt_dXxma = sparse( zeros(nl,nXxma) ); 
    
    for k=1:nXxma
        PertSel =diagmaAux(:,iXxma(k)); %AAB- Selects the column of diagmaAux representing only the active ma for the specified control

        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing Beq in the Perturbed branch (One ma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.*PertSel); 
        
        %Yf and Yt Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert); %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sf and St evaluated in x+pert
        SfPert = diag(Cf*V) * conj(Yf_Pert * V);
        StPert = diag(Ct*V) * conj(Yt_Pert * V);
        
        %Af and At evaluated in x+pert
        AfPert = SfPert .* conj(SfPert);
        AtPert = StPert .* conj(StPert);  
         
        %Partials of Af and At w.r.t. ma Finite differences f'(x) ~~ ( f(x+pert) - f(x) ) / pert 
        num_dAf_dXxma(:, k) = (AfPert - Af )/ pert; %AAB- Final dSf_dma has a size of [nl, nXxma] 
        num_dAt_dXxma(:, k) = (AtPert - At )/ pert; %AAB- Final dSf_dma has a size of [nl, nXxma] 
    end
    
end