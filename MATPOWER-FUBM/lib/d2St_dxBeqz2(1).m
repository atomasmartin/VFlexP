function [H13, H23, H31, H32, H33] = d2St_dxBeqz2(branch, V, mu, vcart)
%D2ST_DXBEQZ2  Computes 2nd derivatives of complex brch power flow w.r.t. BeqzVa, BeqzVm, VaBeqz, VmBeqz, BeqzBeqz.
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 5th argument. So far only Polar has been coded
%
%   [HBzVa, HBzVm, HVaBz, HVmBz, HBzBz] = D2ST_DXBEQZ2(BRANCH, V, MU)
%   [HBzVa, HBzVm, HVaBz, HVmBz, HBzBz] = D2ST_DXBEQZ2(BRANCH, V, MU, 0)
%
%   Returns 5 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HBzVa, HBzVm, HVaBz, HVmBz, HBzBz] = D2ST_DXBEQZ2(BRANCH, V, MU, 1)
%
%   Not Coded Yet
%
%   Examples:
%       il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
%       [HBzVa, HBzVm, HVaBz, HVmBz, HBzBz] = d2St_dxBeqz2(branch(il,:), V, mu, vcart);
%
%       Here the output matrices correspond to:
%           HBzVa = d/dBeqz (dSt_dVa.'   * mu)
%           HBzVm = d/dBeqz (dSt_dVm.'   * mu)
%           HVaBz = d/dVa   (dSt_dBeqz.' * mu)
%           HVmBz = d/dVm   (dSt_dBeqz.' * mu)
%           HBzBz = d/dBeqz (dSt_dBeqz.' * mu)
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

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
%% default input args
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

%% constants
nl = length(mu); %% number of lines
nb = length(V);  %% number of buses
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Beq
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('d2St_dxBeq2: Derivatives of Flow Limit equations w.r.t Beq in cartasian has not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active Beq 
    BeqAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the seclector
    BeqAux1(iBeqz) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqsel = sparse( diag(BeqAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dBeq = sparse( zeros(nl,nBeqz) );
    dYff_dBeq = sparse( zeros(nl,nBeqz) );
    dYft_dBeq = sparse( zeros(nl,nBeqz) );
    dYtf_dBeq = sparse( zeros(nl,nBeqz) ); 
    
    d2St_dBeqzVa = sparse( zeros(nb,   nBeqz) );
    d2St_dBeqzVm = sparse( zeros(nb,   nBeqz) );
    d2St_dBeqz2  = sparse( zeros(nBeqz,nBeqz) );
    
    for k=1:nBeqz 
        for kk=1:nb %dBeqzVx
            %% Second Derivatives
            Beqsel1=diagBeqsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beq
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. Beq
            dYtt_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytt
            dYff_dBeq(:, k) = sparse((  (1j * Beqsel1 )./ ( (k2.*abs(tap)).^2 )  ));  %AAB- Yff does have Beq, this is the derivative
            dYft_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Yft
            dYtf_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytf

            %Partials of Yf, Yt, w.r.t. Beq
            dYf_dBeq = dYff_dBeq(:, k).* Cf + dYft_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq
            dYt_dBeq = dYtf_dBeq(:, k).* Cf + dYtt_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq   

            %2nd Derivatives of Sbus w.r.t. VxBeq
            d2St_dBeqzVa(kk, k) = ((diag(Ct*V)*conj(dYt_dBeq*dVa(:,kk))   +   (Ct*dVa(:,kk)).*conj(dYt_dBeq*V)).')*mu; %AAB- Final dSt_dBeqVa has a size of [nl, nBeqx] 
            d2St_dBeqzVm(kk, k) = ((diag(Ct*V)*conj(dYt_dBeq*dVm(:,kk))   +   (Ct*dVm(:,kk)).*conj(dYt_dBeq*V)).')*mu; %AAB- Final dSt_dBeqVm has a size of [nl, nBeqx] 
        end
        for kk=1:nBeqz
            %% Second Derivatives %The BeqBeq derivative is zero because in the first derivative the beq is eliminated.
            Beqsel2=BeqAux2(:,iBeqz(kk));                         %AAB- From the zero aux matrix we select the column that we will use so there is only 1 element active at the time. It will be zero, but this is how it would have been obtained if it was not zero
            d2Ytf_dBeq2 = Beqsel2;                                %AAB- must be zero
            d2Ytt_dBeq2 = Beqsel2;                                %AAB- must be zero
            d2Yt_dBeq2 = d2Ytf_dBeq2.*Cf + d2Ytt_dBeq2.*Ct;       %AAB- must be zero
            d2St_dBeqz2(kk,k) = (diag(Ct*V)*conj(d2Yt_dBeq2*V)).'*mu;    %AAB- must be zero           
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
H13 = sparse(d2St_dBeqzVa);
H23 = sparse(d2St_dBeqzVm);
H31 = H13.';
H32 = H23.';
H33 = sparse(d2St_dBeqz2);



