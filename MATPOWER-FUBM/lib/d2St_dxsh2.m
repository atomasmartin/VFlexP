function [H15, H25, H35, H45, H51, H52, H53, H54, H55] = d2St_dxsh2(branch, V, mu, vcart)
%D2ST_DXSH2  Computes 2nd derivatives of complex brch power flow "to" w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh.
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 4th argument. So far only Polar has been coded
%
%   [HShVa, HShVm, HShBz, HShBv, HVaSh, HVmSh, HBzSh, HBvSh, HShSh] = D2ST_DXSH2(BRANCH, V, MU)
%   [HShVa, HShVm, HShBz, HShBv, HVaSh, HVmSh, HBzSh, HBvSh, HShSh] = D2ST_DXSH2(BRANCH, V, MU, 0)
%
%   Returns 9 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq of the product of a vector MU with the 1st partial derivatives of 
%   the complex branch power flows.
%
%   [HShVa, HShVm, HShBz, HShBv, HVaSh, HVmSh, HBzSh, HBvSh, HShSh] = D2ST_DXSH2(BRANCH, V, MU, 1)
%
%   Not Coded Yet
%
%   Examples:
%       il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
%       [HShVa, HShVm, HShBz, HShBv, HVaSh, HVmSh, HBzSh, HBvSh, HShSh] = d2St_dxsh2(branch(il,:), V, mu, vcart);
%
%       Here the output matrices correspond to:
%           HShVa = d/dSh   (dSt_dVa.'   * mu)
%           HShVm = d/dSh   (dSt_dVm.'   * mu)
%           HShBz = d/dSh   (dSt_dBeqz.' * mu)
%           HShBv = d/dSh   (dSt_dBeqv.' * mu)
%           HVaSh = d/dVa   (dSt_dSh.'   * mu)
%           HVmSh = d/dVm   (dSt_dSh.'   * mu)
%           HBzSh = d/dBeqz (dSt_dSh.'   * mu)
%           HBvSh = d/dBeqv (dSt_dSh.'   * mu)
%           HShSh = d/dSh   (dSt_dSh.'   * mu)
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
if nargin < 4
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

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift

[stat, Cf, Ct, k2, tap, Ys, Bc] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch

%% Calculation of derivatives
if vcart
    error('d2St_dxsh2: Derivatives of Flow Limit equations w.r.t Theta_sh in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dsh = sparse( zeros(nl,nPfsh) );
    dYff_dsh = sparse( zeros(nl,nPfsh) );
    dYft_dsh = sparse( zeros(nl,nPfsh) );
    dYtf_dsh = sparse( zeros(nl,nPfsh) ); 
    
    d2St_dshVa = sparse( zeros(nb,nPfsh) ); 
    d2St_dshVm = sparse( zeros(nb,nPfsh) );
    d2St_dshBeqz = sparse( zeros(nBeqz,nPfsh) );
    d2St_dshBeqv = sparse( zeros(nBeqv,nPfsh) );
    d2St_dsh2 = sparse( zeros(nPfsh,nPfsh) );
    
    for k=1:nPfsh
        for kk=1:nb %dshVx
            %% Second Derivatives
            Yssh=diagYssh(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Sh
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. Theta_sh
            dYff_dsh(:, k) = sparse( zeros(nl,1) );
            dYft_dsh(:, k) = sparse( -Yssh./(-1j*k2.*conj(tap)) ); %AAB- It also could be: sparse( ( -1j .* Yssh ) ./ ( k2 .* conj(tap) ) );
            dYtf_dsh(:, k) = sparse( -Yssh./( 1j*k2.*tap      ) ); %AAB- It also could be: sparse( (  1j .* Yssh ) ./ ( k2 .*      tap  ) );
            dYtt_dsh(:, k) = sparse( zeros(nl,1) );
        
            %Partials of Yf, Yt, w.r.t. Theta_shift
            %dYf_dsh = dYff_dsh(:, k).* Cf + dYft_dsh(:, k).* Ct; %AAB- size [nl,nb] per active Theta_Shift
            dYt_dsh = dYtf_dsh(:, k).* Cf + dYtt_dsh(:, k).* Ct; %AAB- size [nl,nb] per active Theta_Shift

            %2nd Derivatives of Sbus w.r.t. ShVx
            d2St_dshVa(kk, k) = ((diag(Ct*V)*conj(dYt_dsh*dVa(:,kk))   +   (Ct*dVa(:,kk)).*conj(dYt_dsh*V)).')*mu; %AAB- Final dSf_dShVa has a size of [nl, nPfsh] 
            d2St_dshVm(kk, k) = ((diag(Ct*V)*conj(dYt_dsh*dVm(:,kk))   +   (Ct*dVm(:,kk)).*conj(dYt_dsh*V)).')*mu; %AAB- Final dSf_dShVm has a size of [nl, nPfsh] 
        end
        for kk=1:nBeqz
            %% Second Derivatives %The shBeqz derivative is zero because Yff, Yft, Ytf and Ytt do not share Theta_shift and Beqz for any case.
            %d2Yff_dshBeqz = zeros(nl,1);                                  %AAB- must be zero
            %d2Yft_dshBeqz = zeros(nl,1);                                  %AAB- must be zero
            d2Ytf_dshBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dshBeqz = zeros(nl,1);                                   %AAB- must be zero
            
            %d2Yf_dshBeqz = d2Yff_dshBeqz.*Cf + d2Yft_dshBeqz.*Ct;         %AAB- must be zero
            d2Yt_dshBeqz = d2Ytf_dshBeqz.*Cf + d2Ytt_dshBeqz.*Ct;          %AAB- must be zero
            
            d2St_dshBeqz(kk,k) = (diag(Ct*V)*conj(d2Yt_dshBeqz*V)).'*mu;   %AAB- must be zero  
            
        end
        for kk=1:nBeqv
            %% Second Derivatives %The shBeqv derivative is zero because Yff, Yft, Ytf and Ytt do not share Theta_shift and Beqz for any case.
            %d2Yff_dshBeqv = zeros(nl,1);                                  %AAB- must be zero
            %d2Yft_dshBeqv = zeros(nl,1);                                  %AAB- must be zero
            d2Ytf_dshBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dshBeqv = zeros(nl,1);                                   %AAB- must be zero
            
            %d2Yf_dshBeqv = d2Yff_dshBeqv.*Cf + d2Yft_dshBeqv.*Ct;         %AAB- must be zero
            d2Yt_dshBeqv = d2Ytf_dshBeqv.*Cf + d2Ytt_dshBeqv.*Ct;          %AAB- must be zero
            
            d2St_dshBeqv(kk,k) = (diag(Ct*V)*conj(d2Yt_dshBeqv*V)).'*mu;   %AAB- must be zero  
        end
        for kk=1:nPfsh
            ShSel2=diagShSel(:,iPfsh(kk)); %AAB- Selects the column of diagShSel representing only the active Sh
            
            %% Second Derivatives         
            %d2Yff_dsh2 = sparse( zeros(nl,1) );
            %d2Yft_dsh2 = sparse( (Yssh.*ShSel2)./(k2.*conj(tap)) );        %AAB- Only when Shifter selector 1 and Shifter selector 2 match there will be a derivative (When k = kk). Otherwise zero 
            d2Ytf_dsh2 = sparse( (Yssh.*ShSel2)./(k2.*tap      ) );        %AAB- Only when Shifter selector 1 and Shifter selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Ytt_dsh2 = sparse( zeros(nl,1) );
            
            %d2Yf_dsh2 = d2Yff_dsh2.*Cf + d2Yft_dsh2.*Ct;          
            d2Yt_dsh2 = d2Ytf_dsh2.*Cf + d2Ytt_dsh2.*Ct;          
            
            d2St_dsh2(kk,k) = (diag(Ct*V)*conj(d2Yt_dsh2*V)).'*mu;     %AAB- must be zero           
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
H15 = sparse(d2St_dshVa);
H25 = sparse(d2St_dshVm);
H35 = sparse(d2St_dshBeqz);
H45 = sparse(d2St_dshBeqv);
H51 = H15.';
H52 = H25.';
H53 = H35.';
H54 = H45.';
H55 = sparse(d2St_dsh2);



