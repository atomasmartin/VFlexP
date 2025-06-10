function [H14, H24, H34, H41, H42, H43, H44] = d2Sf_dxBeqv2(branch, V, mu, vcart)
%D2SF_DXBEQV2  Computes 2nd derivatives of complex brch power flow "from" w.r.t. BeqvVa, BeqvVm, BeqvBeqz, VaBeqv, VmBeqv, BeqzBeqv, BeqvBeqv.
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 4th argument. So far only Polar has been coded
%
%   [HBvVa, HBvVm, HBvBz, HVaBv, HVmBv, HBzBv, HBvBv] = D2SF_DXBEQV2(BRANCH, V, MU)
%   [HBvVa, HBvVm, HBvBz, HVaBv, HVmBv, HBzBv, HBvBv] = D2SF_DXBEQV2(BRANCH, V, MU, 0)
%
%   Returns 7 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HBvVa, HBvVm, HBvBz, HVaBv, HVmBv, HBzBv, HBvBv] = D2SF_DXBEQV2(BRANCH, V, MU, 1)
%
%   Not Coded Yet
%
%   Examples:
%       il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
%       [HBvVa, HBvVm, HBvBz, HVaBv, HVmBv, HBzBv, HBvBv] = d2Sf_dxBeqv2(branch(il,:), V, mu, vcart);
%
%       Here the output matrices correspond to:
%           HBvVa = d/dBeqv (dSf_dVa.'   * mu)
%           HBvVm = d/dBeqv (dSf_dVm.'   * mu)
%           HBvBz = d/dBeqv (dSf_dBeqz.' * mu)
%           HVaBv = d/dVa   (dSf_dBeqv.' * mu)
%           HVmBv = d/dVm   (dSf_dBeqv.' * mu)
%           HBzBv = d/dBeqz (dSf_dBeqv.' * mu)
%           HBvBv = d/dBeqv (dSf_dBeqv.' * mu)
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

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC size[nBeqv,1]
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq


[stat, Cf, Ct, k2, tap] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('d2Sf_dxBeq2: Derivatives of Flow Limit equations w.r.t Beq in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active Beq 
    BeqAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the seclector
    BeqAux1(iBeqv) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqsel = sparse( diag(BeqAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dBeq = sparse( zeros(nl,nBeqv) );
    dYff_dBeq = sparse( zeros(nl,nBeqv) );
    dYft_dBeq = sparse( zeros(nl,nBeqv) );
    dYtf_dBeq = sparse( zeros(nl,nBeqv) ); 
    
    d2Sf_dBeqvVa   = sparse( zeros(nb,   nBeqv) ); 
    d2Sf_dBeqvVm   = sparse( zeros(nb,   nBeqv) );
    d2Sf_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    d2Sf_dBeqv2    = sparse( zeros(nBeqv,nBeqv) );
    
    for k=1:nBeqv 
        for kk=1:nb %dBeqzVx
            %% Second Derivatives
            Beqsel1=diagBeqsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. Beqv
            dYtt_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytt
            dYff_dBeq(:, k) = sparse((  (1j * Beqsel1 )./ ( (k2.*abs(tap)).^2 )  ));   %AAB- Yff does have Beq, this is the derivative
            dYft_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Yft
            dYtf_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytf

            %Partials of Yf, Yt, w.r.t. Beqv
            dYf_dBeq = dYff_dBeq(:, k).* Cf + dYft_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq
            dYt_dBeq = dYtf_dBeq(:, k).* Cf + dYtt_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq   

            %2nd Derivatives of Sbus w.r.t. VxBeqv
            d2Sf_dBeqvVa(kk, k) = ((diag(Cf*V)*conj(dYf_dBeq*dVa(:,kk))   +   (Cf*dVa(:,kk)).*conj(dYf_dBeq*V)).')*mu; %AAB- Final dSf_dVaBeq has a size of [nl, nBeqx] 
            d2Sf_dBeqvVm(kk, k) = ((diag(Cf*V)*conj(dYf_dBeq*dVm(:,kk))   +   (Cf*dVm(:,kk)).*conj(dYf_dBeq*V)).')*mu; %AAB- Final dSf_dVmBeq has a size of [nl, nBeqx] 
        end
        for kk=1:nBeqz
            %% Second Derivatives %The BeqvBeqz derivative is zero because a VSC type 1 cannot be a type 2 at the same time.
            d2Yff_dBeqvBeqz = zeros(nl,1);                                    %AAB- must be zero
            d2Yf_dBeqvBeqz = d2Yff_dBeqvBeqz.*Cf;                             %AAB- must be zero
            d2Sf_dBeqvBeqz(kk,k) = (diag(Cf*V)*conj(d2Yf_dBeqvBeqz*V)).'*mu;  %AAB- must be zero  
        end
        
        for kk=1:nBeqv
            %% Second Derivatives %The BeqBeq derivative is zero because in the first derivative the beq is eliminated.
            Beqsel2=BeqAux2(:,iBeqv(kk));                                     %AAB- From the zero aux matrix we select the column that we will use so there is only 1 element active at the time. It will be zero, but this is how it would have been obtained if it was not zero
            d2Yff_dBeq2 = (1j.*Beqsel2)./ ((k2.*abs(tap)).^2);                %AAB- must be zero
            d2Yf_dBeq2 = d2Yff_dBeq2.*Cf;                                     %AAB- must be zero
            d2Sf_dBeqv2(kk,k) = (diag(Cf*V)*conj(d2Yf_dBeq2*V)).'*mu;         %AAB- must be zero           
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
H14 = sparse(d2Sf_dBeqvVa);
H24 = sparse(d2Sf_dBeqvVm);
H34 = sparse(d2Sf_dBeqvBeqz);
H41 = H14.';
H42 = H24.';
H43 = H34.';
H44 = sparse(d2Sf_dBeqv2);



