function [num_dPfdp_dmax] = dPfdp_dmaPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DSBR_DMAPERT   Computes partial derivatives of the branch power flows w.r.t. ma (tap) (Finite differences method).
%
%   ma can be used either to control the Qf, Qt, Vf or Vt to a certain set 
%   value Qfset, Qtset, Vfset or Vtset respectively. Thus, the derivatives 
%   are separated for each function. The derivatives w.r.t. ma for a
%   certain control will be chosen depending on the 3rd argument "ctrl".
%   So that:
%
%   ctrl = 1 : Qf = Qfset, "from side", Transformers (not for VSC)
%   ctrl = 2 : Qt = Qtset, "to side",   Transformers and VSC
%   ctrl = 3 : Vf = Qtset, "from side", Transformers (not for VSC)
%   ctrl = 4 : Vt = Qtset, "to side",   Transformers and VSC
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [NUM_DSF_DQFMAX, NUM_DST_DQFMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 1, PERT, 0)
%   [NUM_DSF_DQFMAX, NUM_DST_DQFMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 1, PERT)
%
%   QF Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Qf control is active for all buses).
%
%   [NUM_DSF_DQTMAX, NUM_DST_DQTMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 2, PERT, 0)
%   [NUM_DSF_DQTMAX, NUM_DST_DQTMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 2, PERT)
%
%   QT Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Qt control is active for all buses).
%
%   [NUM_DSF_DVFMAX, NUM_DST_DVFMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 3, PERT, 0)
%   [NUM_DSF_DVFMAX, NUM_DST_DVFMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 3, PERT)
%
%   VF Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Vf control is active for all buses).
%
%   [NUM_DSF_DVTMAX, NUM_DST_DVTMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 4, PERT, 0)
%   [NUM_DSF_DVTMAX, NUM_DST_DVTMAX] = DSBR_DMAPERT(BASEMVA, BUS, BRANCH, V, 4, PERT)
%
%   VT Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Vt control is active for all buses).
%
%   The following explains the expressions used to form the matrices:
%
%   Sf = diag(Cf*V) * conj(If)         %   St = diag(Ct*V) * conj(It)
%   Sf = diag(Cf*V) * conj(Yf * V)     %   St = diag(Ct*V) * conj(Yt * V)
%
%   where:
%
%       Yf = Yff * Cf + Yft * Ct
%       Yt = Ytf * Cf + Ytt * Ct
%
%       Ytt = Ys + 1j*Bc/2
%       Yff = Gsw+( (Ytt+1j*Beq) ./ ((k2.^2).*tap .* conj(tap))  ) %%<<AAB- FUBM formulation- Original: Yff = Ytt ./ (tap .* conj(tap));
%       Yft = - Ys ./ conj(tap)
%       Ytf = - Ys ./ tap 
%
%       tap = ma .* exp(1j*pi/180 * ShAngle);
%
%   Polar coordinates:
%     Calculation of Ybus, Yf, Yt for Original and perturbed values
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %Original
%       [YbusPert, YfPert, YtPert] = makeYbus(baseMVA, bus, branch_Pert); %Perturbed
%
%     Power Balance Equation Evaluated with original and perturbed values      
%       Sf     = diag(Cf*V) * conj(Yf * V)
%       St     = diag(Ct*V) * conj(Yt * V)
%       SfPert = diag(Cf*V) * conj(YfPert * V)
%       StPert = diag(Ct*V) * conj(YtPert * V)
%
%     Partials of Sbus w.r.t. ma
%       dSf/dma = (SfPert - Sf) / pert
%       dSt/dma = (StPert - St) / pert
%
%   Examples:
%       [num_dSf_dQfmax, num_dSt_dQfmax] = dSbr_dmaPert(baseMVA, bus, branch, V, 1, pert, vcart);
%       [num_dSf_dQtmax, num_dSt_dQtmax] = dSbr_dmaPert(baseMVA, bus, branch, V, 2, pert, vcart);
%       [num_dSf_dVfmax, num_dSt_dVfmax] = dSbr_dmaPert(baseMVA, bus, branch, V, 3, pert, vcart);
%       [num_dSf_dVtmax, num_dSt_dVtmax] = dSbr_dmaPert(baseMVA, bus, branch, V, 4, pert, vcart);
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
%   and Baljinnyam Sereeter, Delft University of Technology
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
    error('dSbr_dmaPert: Control type can only be type 1 (Qf), 2 (Qt), 3(Vf), or 4(Vt)')    
end  

%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nXxma = length(iXxma);      %% AAB- Number of elements with Voltage or Reactive power controlled by ma  
f = branch(:, F_BUS);       %% list of "from" buses
Vm = abs(V);                %% Voltage Magnitude
                              
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dPfdp_dmaPert: Derivatives of Power injections equations w.r.t ma using Finite Differences in cartasian have not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active ma for the specified control 
    maAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    maAux(iXxma) = 1; %AAB- Fill the selector with "1" where ma is active
    diagmaAux = sparse( diag(maAux) );  %AAB- Diagonal of ma selector size [nl,nl]
    
    %Changing Perturbation to Degrees since inside makeYbus it gets changed to radians.
    %pertDeg = (pert*180)/pi;
    
    %Save the Power to control through the branch
    Pfset = branch(:,PF)/baseMVA; %Power constraint for the branch element in p.u.
    
    %Save the Voltage-Droop control settings though the branch (   Pf - Pfset = Kdp*(Vmf - Vmfset)  )
    Kdp    = branch(:,KDP   ); %Voltage Droop Slope   setting for the branch element in p.u.
    Vmfset = branch(:,VF_SET); %Voltage Droop Voltage Setting for the branch element in p.u.
    
    %Yf and Yt Original
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);     %AAB- obtain the Ybus, Yf, Yt
        
    %Sf and St evaluated in x
    Sf = diag(Cf*V) * conj(Yf * V);
    
    %Voltage Droop Equation | Droop Pf(x)   - Pfset = Kdp*(Vmf   - Vmfset)
    Pfdp = -real(Sf)   + Pfset + Kdp.*(Vm(f)  - Vmfset);   %Pfdp(x)
    
    %Dimensionalize (Allocate for computational speed)
    num_dPfdp_dmax = sparse( zeros(nl,nXxma) );  
    
    for k=1:nXxma
        PertSel=diagmaAux(:,iXxma(k)); %AAB- Selects the column of diagshAux representing the location of only the active ma to be perturbed
 
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing ma in the Perturbed branch (One ma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.*PertSel); 
        
        %Yf and Yt Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);     %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sf and St evaluated in x+pert
        SfPert = diag(Cf*V) * conj(Yf_Pert * V);
        
        %Droop evaluated in x+pert
        PfdpPert = -real(SfPert) + Pfset + Kdp.*(Vm(f)  - Vmfset);   %Pfdp(x+h)| Droop Pf(x+h) - Pfset = Kdp*(Vmf   - Vmfset)
        
        %Partials of Pfdp w.r.t. Theta shift Finite differences f'(x) ~~ ( f(x+pert) - f(x) ) / pert         
        num_dPfdp_dmax(:, k) = (PfdpPert - Pfdp )/ pert; %AAB- Final dPfdp_dma has a size of [nl, nXxma]
    end
end