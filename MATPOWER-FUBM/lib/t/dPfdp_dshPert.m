function [num_dPfdp_dshx] = dPfdp_dshPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DPFDP_DSHPERT   Computes partial derivatives of Droop Control w.r.t. Theta_shift. (Finite differences method)
%  
%   Theta_shift can be used either, to control the Pf or Pt to a certain set 
%   value Pfset or Ptset respectively for active power control, or to control
%   the Pf depending on the Voltage Droop Control for VSCIII. Thus, the 
%   derivatives are separated for each function. The derivatives w.r.t. 
%   Theta_shift will be chosen for either Pf, Pt or Droop control, depending 
%   on the 3rd argument "ctrl".
%   So that:
%
%   ctrl = 1 : Pf = Pfset, "from side", Phase Shifter Transformers, VSCI and VSCII
%   ctrl = 2 : Pt = Ptset, "to side", 
%   ctrl = 3 : Pf - Pfset = Kdp*(Vf - Vfset), VSCIII Voltage Droop Control
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   of voltage, depending on the 7th argument. So far only polar
%
%   [NUM_DSF_DSHX, NUM_DST_DSHX] = DSBR_DSHPERT(BASEMVA, BUS, BRANCH, V, 1, PERT, 0)
%
%   PF CONTROL
%   Returns two matrices containing partial derivatives of Sf and St
%   power injections w.r.t Theta_shift,(for all lines).
%
%   [NUM_DSF_DSHX, NUM_DST_DSHX] = DSBR_DSHPERT(BASEMVA, BUS, BRANCH, V, 2, PERT, 0)
%
%   PT CONTROL
%   Returns two matrices containing partial derivatives of Sf and St
%   power injections w.r.t Theta_shift,(for all lines).
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
%     Partials of Sbus w.r.t. shift angle
%       dSf/dsh = (SfPert - Sf) / pert
%       dSt/dsh = (StPert - St) / pert
%
%   Examples:
%       [num_dSf_dshx, num_dSt_dshx] = dSbr_dshPert(baseMVA, bus, branch, V, 1, pert, 0);
%       [num_dSf_dshx, num_dSt_dshx] = dSbr_dshPert(baseMVA, bus, branch, V, 1, pert);
%       [num_dSf_dshx, num_dSt_dshx] = dSbr_dshPert(baseMVA, bus, branch, V, 2, pert, 0);
%       [num_dSf_dshx, num_dSt_dshx] = dSbr_dshPert(baseMVA, bus, branch, V, 2, pert);
%
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

%% selection of Control
if ctrl == 1    %from side
    iPxsh = find( (branch(:,PF  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, SH_MIN )~=-360 | branch(:, SH_MAX )~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4) ); %AAB- Location of the branch elements with Pf control by theta_shift to meet the setting. (Converters and Phase Shifter Transformers, but no VSCIII)
elseif ctrl ==2 %to side
    iPxsh = find( (branch(:,PT  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, SH_MIN )~=-360 | branch(:, SH_MAX )~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4) ); %AAB- Location of the branch elements with Pt control by theta_shift to meet the setting. (Converters and Phase Shifter Transformers, but no VSCIII)
elseif ctrl ==3 %Droop Control
    iPxsh = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
else
    error('dSbr_dshPert: Side for element active power control can only be type 1 (from side), 2 (to side) or 3 (Pf-Vdc Voltage Droop Control)')    
end 
%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nPxsh = length(iPxsh);      %% AAB- Number of elements with active power flow controlled by theta_sh  
f = branch(:, F_BUS);       %% list of "from" buses
Vm = abs(V);                %% Voltage Magnitude
                              
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dPfdp_dshPert: Derivatives of Power injections equations w.r.t Theta_shift using Finite Differences in cartasian have not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active Theta shifters 
    shAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    shAux(iPxsh) = 1; %AAB- Fill the selector with "1" where Theta_shift is active
    diagshAux = sparse( diag(shAux) ); %AAB- Diagonal of Theta_shift selector size [nl,nl]
        
    %Changing Perturbation to Degrees since inside makeYbus it gets changed to radians.
    pertDeg = (pert*180)/pi;
    
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
    num_dPfdp_dshx = sparse( zeros(nl,nPxsh) );  
    
    for k=1:nPxsh
        PertSel=diagshAux(:,iPxsh(k)); %AAB- Selects the column of diagshAux representing the location of only the active shift angles to be perturbed
 
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing Theta_sh in the Perturbed branch (One Shifter at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.*PertSel); 
        
        %Yf and Yt Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);     %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sf and St evaluated in x+pert
        SfPert = diag(Cf*V) * conj(Yf_Pert * V);
        
        %Droop evaluated in x+pert
        PfdpPert_sh = -real(SfPert) + Pfset + Kdp.*(Vm(f)  - Vmfset);   %Pfdp(x+h)| Droop Pf(x+h) - Pfset = Kdp*(Vmf   - Vmfset)
        
        %Partials of Pfdp w.r.t. Theta shift Finite differences f'(x) ~~ ( f(x+pert) - f(x) ) / pert         
        num_dPfdp_dshx(:, k) = (PfdpPert_sh - Pfdp )/ pert; %AAB- Final dPfdp_dsh has a size of [nl, nPxsh]
    end
end