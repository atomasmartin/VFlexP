function [dSf_dshx, dSt_dshx] = dSbr_dsh(branch, V, ctrl, vcart)
%DSBR_DSH   Computes partial derivatives of branch power flows w.r.t. Theta_shift.
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
%   of voltage, depending on the 4th argument. So far only polar
%
%   [DSF_DSHF, DST_DSHF] = DSBR_DSH(BRANCH, V, 1, 0)
%
%   PF CONTROL
%   Returns two matrices containing partial derivatives of Sf and St
%   power injections w.r.t Theta_shift,(for all lines).
%
%   [DSF_DSHT, DST_DSHT] = DSBR_DSH(BRANCH, V, 2, 0)
%
%   PT CONTROL
%   Returns two matrices containing partial derivatives of Sf and St
%   power injections w.r.t Theta_shift,(for all lines).
%
%   [DSF_DSHT, DST_DSHT] = DSBR_DSH(BRANCH, V, 3, 0)
%
%   DROOP CONTROL PF-VF
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
%   Polar coordinates:
%     Partials of Ytt, Yff, Yft and Ytf w.r.t. Theta_shift
%       dYtt/dsh = zeros(nl,1)
%       dYff/dsh = zeros(nl,1)
%       dYft/dsh = -Ys./(-1j*k2.*conj(tap))
%       dYtf/dsh = -Ys./( 1j*k2.*tap      )
%
%     Partials of Yf, Yt, Ybus w.r.t. Theta_shift
%       dYf/dsh = dYff/dsh * Cf + dYft/dsh * Ct
%       dYt/dsh = dYtf/dsh * Cf + dYtt/dsh * Ct
%
%       dYbus/dsh = Cf' * dYf/dsh + Ct' * dYt/dsh    
%
%     Partials of Sbr w.r.t. shift angle
%       dSbr/dsh = diag(V) * conj(dYbus/dsh * V)
%
%   Examples:
%       [dSbr_dshf] = dSbr_dsh(branch, V, 1);
%       [dSbr_dshf] = dSbr_dsh(branch, V, 1, vcart);
%       [dSbr_dsht] = dSbr_dsh(branch, V, 2);
%       [dSbr_dsht] = dSbr_dsh(branch, V, 2, vcart);
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
if nargin < 4
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
    error('dSbr_dsh: Side for element active power control can only be type 1 (from side), 2 (to side) or 3 (Pf-Vdc Voltage Droop Control)')    
end 
%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
nPxsh = length(iPxsh);      %% AAB- Number of elements with active power flow controlled by theta_sh  

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('dSbr_dsh: Derivatives of Power balance equations w.r.t Theta_shift in cartasian have not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active Theta shifters 
    shAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    shAux(iPxsh) = 1; %AAB- Fill the selector with "1" where Theta_shift is active
    diagYssh = sparse( diag(shAux.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dsh = sparse( zeros(nl,nPxsh) );
    dYff_dsh = sparse( zeros(nl,nPxsh) );
    dYft_dsh = sparse( zeros(nl,nPxsh) );
    dYtf_dsh = sparse( zeros(nl,nPxsh) );
    dSf_dshx = sparse( zeros(nl,nPxsh) );
    dSt_dshx = sparse( zeros(nl,nPxsh) );    
    for k=1:nPxsh
        Yssh=diagYssh(:,iPxsh(k)); %AAB- Selects the column of diagYssh representing only the active shift angles
        
        %Partials of Ytt, Yff, Yft and Ytf w.r.t. Theta shift
        dYtt_dsh(:, k) = sparse( zeros(nl,1) );
        dYff_dsh(:, k) = sparse( zeros(nl,1) );
        dYft_dsh(:, k) = sparse( -Yssh./(-1j*k2.*conj(tap)) ); %AAB- It also could be: sparse( ( -1j .* Yssh ) ./ ( k2 .* conj(tap) ) );
        dYtf_dsh(:, k) = sparse( -Yssh./( 1j*k2.*tap      ) ); %AAB- It also could be: sparse( (  1j .* Yssh ) ./ ( k2 .*      tap  ) );

        %Partials of Yf and Yt w.r.t. Theta shift
        dYf_dsh = dYff_dsh(:, k).* Cf + dYft_dsh(:, k).* Ct; %AAB- size [nl,nb] per active Theta shift
        dYt_dsh = dYtf_dsh(:, k).* Cf + dYtt_dsh(:, k).* Ct; %AAB- size [nl,nb] per active Theta shift

        %Partials of Sf and St w.r.t. Theta shift
        dSf_dshx(:, k) = diag(Cf*V) * conj(dYf_dsh * V); %AAB- Final dSf_dsh has a size of [nl, nPxsh]
        dSt_dshx(:, k) = diag(Ct*V) * conj(dYt_dsh * V); %AAB- Final dSt_dsh has a size of [nl, nPxsh]
    end
    
end