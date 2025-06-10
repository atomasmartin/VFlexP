function [dSf_dmax, dSt_dmax] = dSbr_dma(branch, V, ctrl, vcart)
%DSBR_DMA   Computes partial derivatives of the branch power flows w.r.t. ma (tap).
%
%   ma can be used either to control the Qf, Qt, Vf or Vt to a certain set 
%   value Qfset, Qtset, Vfset or Vtset respectively. Thus, the derivatives 
%   are separated for each function. The derivatives w.r.t. ma for a
%   certain control will be chosen depending on the 3rd argument "ctrl".
%   So that:
%
%   ctrl = 1 : Qf = Qfset, "from side", Transformers (not for VSC)
%   ctrl = 2 : Qt = Qtset, "to side",   Transformers and VSC
%   ctrl = 3 : Vf = Vfset, "from side", Transformers (not for VSC)
%   ctrl = 4 : Vt = Vtset, "to side",   Transformers and VSC
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument. So far only polar
%
%   [DSF_DQFMAX, DST_DQFMA] = DSBR_DMA(BRANCH, V, 1)
%   [DSF_DQFMAX, DST_DQFMA] = DSBR_DMA(BRANCH, V, 1, VCART)
%
%   QF Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Qf control is active for all buses).
%
%   [DSF_DQTMAX, DST_DQTMA] = DSBR_DMA(BRANCH, V, 2)
%   [DSF_DQTMAX, DST_DQTMA] = DSBR_DMA(BRANCH, V, 2, VCART)
%
%   QT Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Qt control is active for all buses).
%
%   [DSF_DVFMAX, DST_DVFMA] = DSBR_DMA(BRANCH, V, 3)
%   [DSF_DVFMAX, DST_DVFMA] = DSBR_DMA(BRANCH, V, 3, VCART)
%
%   VF Control
%   Returns one matrix containing partial derivatives of the complex bus
%   power injections w.r.t ma (where Vf control is active for all buses).
%
%   [DSF_DVTMAX, DST_DVTMA] = DSBR_DMA(BRANCH, V, 4)
%   [DSF_DVTMAX, DST_DVTMA] = DSBR_DMA(BRANCH, V, 4, VCART)
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
%     Partials of Ytt, Yff, Yft and Ytf w.r.t. ma
%       dYtt/dma = zeros(nl,1)
%       dYff/dma = -2*Yttma./( (k2.^2).*((abs(tap)).^3) )
%       dYft/dma = Ys./( k2.*(abs(tap).*conj(tap)) )
%       dYtf/dma = Ys./( k2.*(abs(tap).*     tap ) )
%
%     Partials of Yf, Yt, Ybus w.r.t. ma
%       dYf/dma = dYff/dma * Cf + dYft/dma * Ct
%       dYt/dma = dYtf/dma * Cf + dYtt/dma * Ct 
%
%     Partials of Sf and St w.r.t. ma
%       dSf/dma = diag(Cf*V) * conj(dYf/dma * V)
%       dSt/dma = diag(Cf*V) * conj(dYt/dma * V)
%
%   Examples:
%       [dSf_dQfmax, dSt_dQfmax] = dSbr_dma(branch, V, 1, vcart);
%       [dSf_dQtmax, dSt_dQtmax] = dSbr_dma(branch, V, 2, vcart);
%       [dSf_dVfmax, dSt_dVfmax] = dSbr_dma(branch, V, 3, vcart);
%       [dSf_dVtmax, dSt_dVtmax] = dSbr_dma(branch, V, 4, vcart);
%
%       [dSf_dQfmax, dSt_dQfmax] = dSbr_dma(branch, V, 1);
%       [dSf_dQtmax, dSt_dQtmax] = dSbr_dma(branch, V, 2);
%       [dSf_dVfmax, dSt_dVfmax] = dSbr_dma(branch, V, 3);
%       [dSf_dVtmax, dSt_dVtmax] = dSbr_dma(branch, V, 4);
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
YttBeq = Ys + 1j*Bc/2 + 1j*Beq; %Ytt + jBeq
if vcart
    error('dSbr_dma: Derivatives of Power balance equations w.r.t ma in cartasian have not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active ma for the specified control 
    maAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    maAux(iXxma) = 1; %AAB- Fill the selector with "1" where ma is active
    diagYsma = sparse( diag(maAux.*Ys) );  %AAB- ma selector multilied by the series addmitance Ys,  size [nl,nl]
    diagYttBeqma= sparse( diag(maAux.*YttBeq) ); %AAB- ma selector multilied by the series addmitance Ytt, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dma = sparse( zeros(nl,nXxma) );
    dYff_dma = sparse( zeros(nl,nXxma) );
    dYft_dma = sparse( zeros(nl,nXxma) );
    dYtf_dma = sparse( zeros(nl,nXxma) );
    dSf_dmax = sparse( zeros(nl,nXxma) );
    dSt_dmax = sparse( zeros(nl,nXxma) );    
    for k=1:nXxma
        Ysma =diagYsma(:,iXxma(k)); %AAB- Selects the column of diagYsma representing only the active ma for the specified control multiplied by Ys
        YttBeqma=diagYttBeqma(:,iXxma(k)); %AAB- Selects the column of diagmaAux representing only the active ma for the specified control
        
        %Partials of Ytt, Yff, Yft and Ytf w.r.t. ma
        dYtt_dma(:, k) = sparse( zeros(nl,1) );
        dYff_dma(:, k) = sparse( -2*YttBeqma./( (k2.^2).*((abs(tap)).^3) ) );
        dYft_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*conj(tap)) ) ); 
        dYtf_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*     tap ) ) ); 

        %Partials of Yf and Yt w.r.t. ma
        dYf_dma = dYff_dma(:, k).* Cf + dYft_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma
        dYt_dma = dYtf_dma(:, k).* Cf + dYtt_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma

        %Partials of Sf and St w.r.t. ma
        dSf_dmax(:, k) = diag(Cf*V) * conj(dYf_dma * V); %AAB- Final dSf_dma has a size of [nl, nXxma] 
        dSt_dmax(:, k) = diag(Ct*V) * conj(dYt_dma * V); %AAB- Final dSf_dma has a size of [nl, nXxma] 
    end
    
end