function [num_dPfdp_dBeqx] = dPfdp_dBeqPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DPFDP_DBEQPERT   Computes partial derivatives of branch power flows w.r.t. Beq (Finite differences method).
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
%   of voltage, depending on the 4th argument. So far only polar
%
%   [NUM_DSF_DBEQX, NUM_DST_DBEQX] = DSBR_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT, 0)
%   [NUM_DSF_DBEQX, NUM_DST_DBEQX] = DSBR_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT)
%
%   Returns one matrix containing partial derivatives of the branch
%   power injections Sf, St w.r.t Beq, (for all lines).
%
%   [NUM_DSF_DBEQX, NUM_DST_DBEQX] = DSBR_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT, 1)
%
%   Not Coded yet
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
%     Partials of Sbus w.r.t. Beq
%       dSf/dma = (SfPert - Sf) / pert
%       dSt/dma = (StPert - St) / pert
%
%   Examples:
%       [num_dSf_dBeqz, num_dSt_dBeqz] = dSbr_dBeqPert(baseMVA, bus, branch, V, 1, pert, vcart);
%       [num_dSf_dBeqv, num_dSt_dBeqv] = dSbr_dBeqPert(baseMVA, bus, branch, V, 2, pert, vcart);
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
f = branch(:, F_BUS);       %% list of "from" buses
Vm = abs(V);                %% Voltage Magnitude
                              
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dPfdp_dBeqPert: Derivatives of Droop equations w.r.t ma using Finite Differences in cartasian have not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active Beq 
    BeqAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    BeqAux(iBeqx) = 1; %AAB- Fill the selector with 1 where Beq is active
    diagBeqAux = sparse( diag(BeqAux) ); %AAB- Beq Selector [nl,nl]
 
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
    num_dPfdp_dBeqx = sparse( zeros(nl,nBeqx) );  
    
    for k=1:nBeqx
        PertSel=diagBeqAux(:,iBeqx(k)); %AAB- Selects the column of diagshAux representing the location of only the active ma to be perturbed
 
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.*PertSel); 
        
        %Yf and Yt Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);     %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sf and St evaluated in x+pert
        SfPert = diag(Cf*V) * conj(Yf_Pert * V);
        
        %Droop evaluated in x+pert
        PfdpPert = -real(SfPert) + Pfset + Kdp.*(Vm(f)  - Vmfset);   %Pfdp(x+h) | Droop Pf(x+h) - Pfset = Kdp*(Vmf   - Vmfset)
        
        %Partials of Pfdp w.r.t. Theta shift Finite differences f'(x) ~~ ( f(x+pert) - f(x) ) / pert         
        num_dPfdp_dBeqx(:, k) = (PfdpPert - Pfdp )/ pert; %AAB- Final dPfdp_dma has a size of [nl, nXxma]
    end
end