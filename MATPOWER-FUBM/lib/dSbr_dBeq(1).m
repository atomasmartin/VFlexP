function [dSf_dBeqx, dSt_dBeqx] = dSbr_dBeq(branch, V, ctrl, vcart)
%DSBBR_DBEQ   Computes partial derivatives of branch power flows w.r.t. Beq.
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
%   [DSF_DBEQX, DST_DBEQX] = DSBR_DBEQ(BRANCH, V, VSC)
%   [DSF_DBEQX, DST_DBEQX] = DSBR_DBEQ(BRANCH, V, VSC, 0)
%
%   Returns one matrix containing partial derivatives of the branch
%   power injections Sf, St w.r.t Beq, (for all lines).
%
%   [DSF_DBEQX, DST_DBEQX] = DSBR_DBEQ(BRANCH, V, VSC, 1)
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
%   Polar coordinates:
%     Partials of Ytt, Yff, Yft and Ytf w.r.t. Beq
%       dYtt/dBeq = zeros(nl,1)
%       dYff/dBeq = ( (1j * ones(nl,1) ) ./ ((k2.^2).*tap .* conj(tap))  )
%       dYft/dBeq = zeros(nl,1)
%       dYtf/dBeq = zeros(nl,1)
%
%     Partials of Yf, Yt, Ybus w.r.t. Beq
%       dYf/dBeq = dYff/dBeq * Cf + dYft/dBeq * Ct
%       dYt/dBeq = dYtf/dBeq * Cf + dYtt/dBeq * Ct   
%
%     Partials of S w.r.t. Beq
%       dSf/dBeq = diag(Cf*V) * conj(dYf/dBeq * V)
%       dSt/dBeq = diag(Ct*V) * conj(dYt/dBeq * V)
%
%   Examples:
%       [dSf_dBeqx, dSt_dBeqx] = dSbr_dBeq(branch, V, vsc);
%       [dSf_dBeqx, dSt_dBeqx] = dSbr_dBeq(branch, V, vsc, vcart);
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
    error('dSbr_dBeq: Derivatives of Power balance equations w.r.t Beq in cartasian has not been coded yet')    

else %AAB- Polar Version
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    
    %Selector of active Beq 
    BeqAux = zeros(nl,1);%AAB- Vector of zeros for the seclector
    BeqAux(iBeqx) = 1; %AAB- Fill the selector with 1 where Beq is active
    diagBeqsel = sparse( diag(BeqAux) ); %AAB- Beq Selector [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dBeq = sparse( zeros(nl,nBeqx) );
    dYff_dBeq = sparse( zeros(nl,nBeqx) );
    dYft_dBeq = sparse( zeros(nl,nBeqx) );
    dYtf_dBeq = sparse( zeros(nl,nBeqx) );
    dSf_dBeqx = sparse( zeros(nl,nBeqx) );
    dSt_dBeqx = sparse( zeros(nl,nBeqx) );
    
    for k=1:nBeqx
        Beqsel=diagBeqsel(:,iBeqx(k)); %AAB- Selects the column of diagBeqsel representing only the active Beq
        
        %Partials of Ytt, Yff, Yft and Ytf w.r.t. Beq
        dYtt_dBeq(:, k) = sparse( zeros(nl,1) );
        dYff_dBeq(:, k) = sparse((  (1j * Beqsel )./ ( (k2.*abs(tap)).^2 )  ));
        dYft_dBeq(:, k) = sparse( zeros(nl,1) );
        dYtf_dBeq(:, k) = sparse( zeros(nl,1) );

        %Partials of Yf, Yt, Ybus w.r.t. Beq
        dYf_dBeq = dYff_dBeq(:, k).* Cf + dYft_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq
        dYt_dBeq = dYtf_dBeq(:, k).* Cf + dYtt_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq
        
        %Partials of Sf and St w.r.t. Beq
        dSf_dBeqx(:, k) = diag(Cf*V) * conj(dYf_dBeq * V); %AAB- Final dSf_dBeq has a size of [nl, nBeqx]
        dSt_dBeqx(:, k) = diag(Ct*V) * conj(dYt_dBeq * V); %AAB- Final dSt_dBeq has a size of [nl, nBeqx]
    end
    
end