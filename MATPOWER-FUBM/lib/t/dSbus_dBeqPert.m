function [num_dSbus_dBeqx] = dSbus_dBeqPert(baseMVA, bus, branch, V, ctrl, pert, vcart)
%DSBUS_DBEQPERT   Computes partial derivatives of power injection w.r.t. Beq.   (Finite differences method)
%
%   Beq can be used either to control the Vdc to a certain set value Vfset 
%   or the Qf to match zero (zero constraint). So the derivatives are
%   separated for each function. The derivatives w.r.t. Beq will be 
%   chosen for either VSC type I and IIIz, or type II, or VSC type I and IIIz 
%   and III depending on the 3rd argument.
%   So that:
%
%   ctrl = 1 : Qf = 0,    Zero constraint  VSCI and VSCIIIz only for Power flow
%   ctrl = 2 : Vf = Vset, Vdc control      VSCII
%   ctrl = 3 : Qf = 0,    Zero constraint  VSCI, VSCIIIz, VSCIII only for OPF
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument. So far only polar
%
%   [NUM_DSBUS_DBEQX] = DSBUS_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT, 0)
%   [NUM_DSBUS_DBEQX] = DSBUS_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT)
%
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t Beq, respectively (for all buses).
%
%   [NUM_DSBUS_DBEQX] = DSBUS_DBEQPERT(BASEMVA, BUS, BRANCH, V, VSC, PERT, 1)
%
%   Not Coded yet
%
%   The following explains the expressions used to form the matrices:
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%   S = diag(V) * conj(Ybus * V)
%   where:
%       Ybus = Cf' * Yf + Ct' * Yt + diag(Ysh)
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
%     Calculation of Ybus for Original and perturbed values
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %Original
%       [YbusPert, YfPert, YtPert] = makeYbus(baseMVA, bus, branch_Pert); %Perturbed
%
%     Power Balance Equation Evaluated with original and perturbed values      
%       Sbus     = diag(V) * conj(Ybus     * V)
%       SbusPert = diag(V) * conj(YbusPert * V)
%
%     Partials of Sbus w.r.t. Beq Finite Differences Method f'(x) ~~ ( f(x+pert) - f(x) ) / pert 
%       dSbus/dsh = (Sbus - SbusPert) / pert
%
%   Examples:
%       [num_dSbus_dBeqz] = dSbus_dBeqPert(baseMVA, bus, branch, V, 1, pert, 0);
%       [num_dSbus_dBeqv] = dSbus_dBeqPert(baseMVA, bus, branch, V, 2, pert, 0);
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
if ctrl == 1 %VSC I and VSCIIIz POWER FLOW
    iBeqx = find ( ( branch(:,CONV) == 1 | branch(:,CONV) == 3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC size[nBeqz,1]
elseif ctrl ==2 %VSC II
    iBeqx = find (branch(:,CONV) == ctrl & branch(:, BR_STATUS)==1 &  branch(:, VF_SET)~=0) ; %AAB- Find branch locations of VSC size[nBeqv,1]
elseif ctrl == 3 %VSC I, VSCIIIz, VSCIII, OPTIMAL POWER FLOW
    iBeqx = find ( ( branch(:,CONV) == 1 | branch(:,CONV) == 3 | branch(:,CONV) == 4 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC size[nBeqz,1]
else
    error('dSbus_dBeq: VSC can only be control 1 (VSCI and VSCIIIz), 2 (VSCII), OR 3 (VSCI, VSCIIIz and VSCIII)')    
end  
%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nBeqx = length(iBeqx);      %% AAB- Number of VSC with active Beq

if vcart
    error('dSbus_dBeqPert: Derivatives of Power balance equations w.r.t Beq using Finite Differences in cartasian has not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active Beq 
    BeqAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    BeqAux(iBeqx) = 1; %AAB- Fill the selector with 1 where Beq is active
    diagBeqAux = sparse( diag(BeqAux) ); %AAB- Beq Selector [nl,nl]
    
    %Ybus Original
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);     %AAB- obtain the Perturbed Ybus and Y from for the first element.
        
    %Sbus evaluated in x
    Sbus = diagV * conj(Ybus * V);    
    
    %Dimensionalize (Allocate for computational speed)
    num_dSbus_dBeqx = sparse( zeros(nb,nBeqx) );
    
    for k=1:nBeqx
        PertSel = diagBeqAux(:,iBeqx(k)); %AAB- Selects the column of diagBeqsel representing only the active Beq

        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 

        %Ybus Perturbed
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);     %AAB- obtain the Perturbed Ybus, Yf, and Yt from for the k element.
        
        %Sbus evaluated in x+pert
        SbusPert = diagV * conj(Ybus_Pert * V);
        
        %Partials of S w.r.t. Beq
        num_dSbus_dBeqx(:, k) = (SbusPert - Sbus )/ pert; %AAB- Final dSbus_dBeq has a size of [nb, nBeqx] 
    end
end