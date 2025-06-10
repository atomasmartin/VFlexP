function [G16, G26, G56, G61, G62, G65, G66] = d2Sbus_dxBeqv2(branch, V, lam, vcart)
%D2SBUS_DXBEQV2   Computes 2nd derivatives of power injection w.r.t. BeqvVa, BeqvVm, BeqvBz, VaBeqv, VmBeqv, BeqzBeqv, BeqvBeqv.
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument. So far only polar
%
%   [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = D2SBUS_DXBEQV2(BRANCH, V, LAM)
%   [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = D2SBUS_DXBEQV2(BRANCH, V, LAM, 0)
%
%   Returns 7 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. BeqvVa, BeqvVm, BeqvBz, VaBeqv, VmBeqv, BeqzBeqv,
%   and BeqvBeqv
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GBvVa, GBvVm, GBvBz, GVaBv, GVmBv, GBzBv, GBvBv] = d2Sbus_dxBeqv2(branch, V, lam, vcart);
%
%       Here the output matrices correspond to:
%           GBvVa = d/dVa   (dSbus_dBeqv.' * lam)
%           GBvVm = d/dVm   (dSbus_dBeqv.' * lam)
%           GBvBz = d/dBeqz (dSbus_dBeqv.' * lam)
%           GVaBv = d/dBeqv (dSbus_dVa.'   * lam)
%           GVmBv = d/dBeqv (dSbus_dVm.'   * lam)
%           GBzBv = d/dBeqv (dSbus_dBeqz.' * lam)
%           GBvBv = d/dBeqv (dSbus_dBeqv.' * lam)

%
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
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
if nargin < 4
    vcart = 0;      %% default to polar coordinates
end

%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines
%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC size[nBeqv,1]
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

if vcart
    error('d2Sbus_dxBeq2: Derivatives of Power balance equations w.r.t Beq in cartasian has not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active Beqv 
    BeqAux1 = sparse( zeros(nl,1) );      %AAB- Vector of zeros for the selector
    BeqAux1(iBeqv) = 1;                   %AAB- Fill the selector with 1 where Beq is active
    diagBeqsel = sparse( diag(BeqAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqAux2 = sparse( zeros(nl,nl));      %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dBeq = sparse( zeros(nl,nBeqv) );
    dYff_dBeq = sparse( zeros(nl,nBeqv) );
    dYft_dBeq = sparse( zeros(nl,nBeqv) );
    dYtf_dBeq = sparse( zeros(nl,nBeqv) ); 
    
    d2Sbus_dBeqvVa = sparse( zeros(nb,nBeqv) ); 
    d2Sbus_dBeqvVm = sparse( zeros(nb,nBeqv) ); 
    d2Sbus_dBeqvBeqz = sparse( zeros(nBeqz,nBeqv) );
    d2Sbus_dBeqv2 = sparse( zeros(nBeqv,nBeqv) );
    
    for k=1:nBeqv 
        for kk=1:nb %dBeqvVx
            %% Second Derivatives    
            Beqsel1=diagBeqsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. Beqv
            dYtt_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytt
            dYff_dBeq(:, k) = sparse((  (1j * Beqsel1 )./ ( (k2.*abs(tap)).^2 )  ));   %AAB- Yff does have Beq, this is the derivative
            dYft_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Yft
            dYtf_dBeq(:, k) = sparse( zeros(nl,1) );                                  %AAB- Only zeros because there is no Beq in Ytf

            %Partials of Yf, Yt, Ybus w.r.t. Beqv
            dYf_dBeq = dYff_dBeq(:, k).* Cf + dYft_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq
            dYt_dBeq = dYtf_dBeq(:, k).* Cf + dYtt_dBeq(:, k).* Ct; %AAB- size [nl,nb] per active Beq

            dYbus_dBeq = Cf' * dYf_dBeq + Ct' * dYt_dBeq;     %AAB- size [nb,nb] per active Beq       

            %2nd Derivatives of Sbus w.r.t. VxBeqv
            d2Sbus_dBeqvVa(kk, k) = ((dVa(:,kk).*conj(dYbus_dBeq*V) + V.*conj(dYbus_dBeq*(dVa(:,kk)))).')*lam;       %AAB- Final d2Sbus_dBeqVa has a size of [nb, nBeqx] 
            d2Sbus_dBeqvVm(kk, k) = ((dVm(:,kk).*conj(dYbus_dBeq*V) + V.*conj(dYbus_dBeq*(dVm(:,kk)))).')*lam; %AAB- Final d2Sbus_dBeqVm has a size of [nb, nBeqx] 

        end
        
        for kk=1:nBeqz
            %% Second Derivatives %The BeqvBeqz derivative is zero because a VSC type 1 cannot be a type 2 at the same time.
            d2Yff_dBeqvBeqz = zeros(nl,1);                                 %AAB- must be zero
            d2Ybus_dBeqvBeqz = Cf'*(d2Yff_dBeqvBeqz.*Cf);                  %AAB- must be zero
            d2Sbus_dBeqvBeqz(kk,k) = (V.*conj(d2Ybus_dBeqvBeqz*V)).'*lam;  %AAB- must be zero           
        end
        for kk=1:nBeqv
            %% Second Derivatives %The BeqvBeqv derivative is zero because in the first derivative the beq is eliminated.
            Beqsel2=BeqAux2(:,iBeqv(kk));                                  %AAB- From the zero aux matrix we select the column that we will use so there is only 1 element active at the time. It will be zero, but this is how it would have been obtained if it was not zero
            d2Yff_dBeq2 = (1j.*Beqsel2)./ ((k2.*abs(tap)).^2);             %AAB- must be zero
            d2Ybus_dBeq2 = Cf'*(d2Yff_dBeq2.*Cf);                          %AAB- must be zero
            d2Sbus_dBeqv2(kk,k) = (V.*conj(d2Ybus_dBeq2*V)).'*lam;         %AAB- must be zero    
            
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
G16 = sparse(d2Sbus_dBeqvVa);
G26 = sparse(d2Sbus_dBeqvVm);
G56 = sparse(d2Sbus_dBeqvBeqz);
G61 = G16.';
G62 = G26.';
G65 = G56.';
G66 = sparse(d2Sbus_dBeqv2);



