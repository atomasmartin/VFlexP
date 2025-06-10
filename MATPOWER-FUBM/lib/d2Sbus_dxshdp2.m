function [G110, G210, G510, G610, G710, G810, G910, G101, G102, G105, G106, G107, G108, G109, G1010] = d2Sbus_dxshdp2(branch, V, lam, vcart)
%D2SBUS_DXSHDP2   Computes 2nd derivatives of power injection w.r.t. ShdpVa, ShdpVm, ShdpBeqz, ShdpBeqv, ShdpSh, Shdpqtma, VaShdp, VmShdp, BeqzShdp, BeqvShdp, ShShdp, qtmaShdp, vtmaShdp, ShdpShdp
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument. So far only polar
%
%   [GShdpVa, GShdpVm, GShdpBeqz, GShdpBeqv, GShdpSh, GShdpqtma, GShdpvtma, GVaShdp, GVmShdp, GBeqzShdp, GBeqvShdp, GShShdp, GqtmaShdp, GvtmaShdp, GShdpShdp] = D2SBUS_DXSHDP2(BRANCH, V, LAM)

%   [GShdpVa, GShdpVm, GShdpBeqz, GShdpBeqv, GShdpSh, GShdpqtma, GShdpvtma, GVaShdp, GVmShdp, GBeqzShdp, GBeqvShdp, GShShdp, GqtmaShdp, GvtmaShdp, GShdpShdp] = D2SBUS_DXSHDP2(BRANCH, V, LAM, 0)
%
%   Returns 15 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. ShdpVa, ShdpVm, ShdpBeqz, ShdpBeqv, ShdpSh, Shdpqtma, Shdpvtma, 
%   VaShdp, VmShdp, BeqzShdp, BeqvShdp, ShShdp, qtmaShdp, vtmaShdp and ShdpShdp.
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GShdpVa, GShdpVm, GShdpBeqz, GShdpBeqv, GShdpSh, GShdpqtma, GShdpvtma,...
%                 GVaShdp, GVmShdp, GBeqzShdp, GBeqvShdp, GShShdp, GqtmaShdp, GvtmaShdp, GShdpShdp] = ...
%                                                             d2Sbus_dxqtma2(branch, V, lam, vcart);
%
%       Here the output matrices correspond to:
%           GShdpVa   = d/dVa   (dSbus_dShdp.' * lam)
%           GShdpVm   = d/dVm   (dSbus_dShdp.' * lam)
%           GShdpBz   = d/dBeqz (dSbus_dShdp.' * lam)
%           GShdpBv   = d/dBeqv (dSbus_dShdp.' * lam)
%           GShdpSh   = d/dsh   (dSbus_dShdp.' * lam)
%           GShdpqtma = d/dqtma (dSbus_dShdp.' * lam)
%           GShdpvtma = d/dvtma (dSbus_dShdp.' * lam)
%           GVaShdp   = d/dShdp (dSbus_dVa.'   * lam)
%           GVmShdp   = d/dShdp (dSbus_dVm.'   * lam)
%           GBzShdp   = d/dShdp (dSbus_dBeqz.' * lam)
%           GBvShdp   = d/dShdp (dSbus_dBeqv.' * lam)
%           GShShdp   = d/dShdp (dSbus_dsh.'   * lam)
%           GqtmaShdp = d/dShdp (dSbus_dqtma.' * lam)
%           GvtmaShdp = d/dShdp (dSbus_dvtma.' * lam)
%           GShdpShdp = d/dShdp (dSbus_dShdp.' * lam)
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
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxshdp2: Derivatives of Power balance equations w.r.t ma/tap in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active Theta_dp Droop 
    QtmaSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector ma controlling Qt
    VtmaSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector ma controlling Vt
    BeqzSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqz controlling Zero Constraint
    BeqvSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqv controlling Vf
    PfshSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector sh controlling Pf  
    PfdpSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector shdp controlling Pf  
    
    QtmaSel(iQtma) = 1;                        %AAB- Fill the selector with 1 where ma is active and controlling Qt
    VtmaSel(iVtma) = 1;                        %AAB- Fill the selector with 1 where ma is active and controlling Vt
    BeqzSel(iBeqz) = 1;                        %AAB- Fill the selector with 1 where Beqz is active and controlling Zero Constraint
    BeqvSel(iBeqv) = 1;                        %AAB- Fill the selector with 1 where Beqv is active and controlling Vf
    PfshSel(iPfsh) = 1;                        %AAB- Fill the selector with 1 where sh is active and controlling Pf
    PfdpSel(iPfdp) = 1;                        %AAB- Fill the selector with 1 where shdp is active and controlling Pf
    
    diagQtmaSel = sparse( diag(QtmaSel));      %AAB- Diagonal of the selector where ma is active and controlling Qt
    diagVtmaSel = sparse( diag(VtmaSel));      %AAB- Diagonal of the selector where ma is active and controlling Qt
    diagPfshSel = sparse( diag(PfshSel));      %AAB- Diagonal of the selector where sh is active and controlling Pf
    diagBeqzSel = sparse( diag(BeqzSel));      %AAB- Diagonal of the selector where Beqz is active and controlling Zero Constraint
    diagBeqvSel = sparse( diag(BeqvSel));      %AAB- Diagonal of the selector where Beqv is active and controlling Vf   
    diagPfdpSel = sparse( diag(PfdpSel));      %AAB- Diagonal of the selector where shdp is active and controlling Pf
    
    diagYsshdp   = sparse( diag(PfdpSel.*Ys) );        %AAB- Theta_Droop selector multilied by the series addmitance Ys, size [nl,nl]
    diagYsVtma   = sparse( diag(VtmaSel.*Ys) );        %AAB- ma/tap selector multilied by the series addmitance Ys, size [nl,nl]
    diagYttBVtma = sparse( diag(VtmaSel.*YttB) );      %AAB- ma selector multilied by the series addmitance Ytt, size [nl,nl]

    %Dimensionalize (Allocate for computational speed)
    dYtt_dshdp = sparse( zeros(nl,nPfdp) );
    dYff_dshdp = sparse( zeros(nl,nPfdp) );
    dYft_dshdp = sparse( zeros(nl,nPfdp) );
    dYtf_dshdp = sparse( zeros(nl,nPfdp) ); 
    
    d2Sbus_dshdpVa   = sparse( zeros(nb,   nPfdp) ); 
    d2Sbus_dshdpVm   = sparse( zeros(nb,   nPfdp) ); 
    d2Sbus_dshdpBeqz = sparse( zeros(nBeqz,nPfdp) );
    d2Sbus_dshdpBeqv = sparse( zeros(nBeqv,nPfdp) );
    d2Sbus_dshdpsh   = sparse( zeros(nPfsh,nPfdp) );
    d2Sbus_dshdpqtma = sparse( zeros(nQtma,nPfdp) );
    d2Sbus_dshdpvtma = sparse( zeros(nVtma,nPfdp) );    
    d2Sbus_dshdp2    = sparse( zeros(nPfdp,nPfdp) );
    
    for k=1:nPfdp
        for kk=1:nb %dPfdpVx
            %% Second Derivatives  
            Ysshdp=diagYsshdp(:,iPfdp(k)); %AAB- Selects the column of diagShsel representing only the active Sh
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. Theta_dp
            dYff_dshdp(:, k) = sparse( zeros(nl,1) );
            dYft_dshdp(:, k) = sparse( -Ysshdp./(-1j*k2.*conj(tap)) ); %AAB- It also could be: sparse( ( -1j .* Yssh ) ./ ( k2 .* conj(tap) ) );
            dYtf_dshdp(:, k) = sparse( -Ysshdp./( 1j*k2.*tap      ) ); %AAB- It also could be: sparse( (  1j .* Yssh ) ./ ( k2 .*      tap  ) );
            dYtt_dshdp(:, k) = sparse( zeros(nl,1) );
        
            %Partials of Yf, Yt, w.r.t. Theta_shift Droop
            dYf_dshdp = dYff_dshdp(:, k).* Cf + dYft_dshdp(:, k).* Ct; %AAB- size [nl,nb] per active Theta_Shift
            dYt_dshdp = dYtf_dshdp(:, k).* Cf + dYtt_dshdp(:, k).* Ct; %AAB- size [nl,nb] per active Theta_Shift

            dYbus_dshdp = Cf' * dYf_dshdp + Ct' * dYt_dshdp;     %AAB- size [nb,nb] per active ma       

            %2nd Derivatives of Sbus w.r.t. shdpVx
            d2Sbus_dshdpVa(kk, k) = ((dVa(:,kk).*conj(dYbus_dshdp*V) + V.*conj(dYbus_dshdp*(dVa(:,kk)))).')*lam; %AAB- Final d2Sbus_dshdpVa has a size of [nb, nshdp] 
            d2Sbus_dshdpVm(kk, k) = ((dVm(:,kk).*conj(dYbus_dshdp*V) + V.*conj(dYbus_dshdp*(dVm(:,kk)))).')*lam; %AAB- Final d2Sbus_dshdpVm has a size of [nb, nshdp] 
        end
        for kk=1:nBeqz
            %% Second Derivatives %The shBeqz derivative is zero because Yff, Yft, Ytf and Ytt do not share Theta_shift and Beqz for any case.
            d2Yff_dshdpBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Yft_dshdpBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytf_dshdpBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dshdpBeqz = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dshdpBeqz = d2Yff_dshdpBeqz.*Cf + d2Yft_dshdpBeqz.*Ct;      %AAB- must be zero
            d2Yt_dshdpBeqz = d2Ytf_dshdpBeqz.*Cf + d2Ytt_dshdpBeqz.*Ct;      %AAB- must be zero
            
            d2Ybus_dshdpBeqz = Cf'*(d2Yf_dshdpBeqz) + Ct'*(d2Yt_dshdpBeqz);  %AAB- must be zero
            d2Sbus_dshdpBeqz(kk,k) = (V.*conj(d2Ybus_dshdpBeqz*V)).'*lam;    %AAB- must be zero           
        end
        for kk=1:nBeqv
            %% Second Derivatives %The shBeqv derivative is zero because Yff, Yft, Ytf and Ytt do not share Theta_shift and Beqv for any case.
            d2Yff_dshdpBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Yft_dshdpBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytf_dshdpBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dshdpBeqv = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dshdpBeqv = d2Yff_dshdpBeqv.*Cf + d2Yft_dshdpBeqv.*Ct;      %AAB- must be zero
            d2Yt_dshdpBeqv = d2Ytf_dshdpBeqv.*Cf + d2Ytt_dshdpBeqv.*Ct;      %AAB- must be zero
            
            d2Ybus_dshdpBeqv = Cf'*(d2Yf_dshdpBeqv) + Ct'*(d2Yt_dshdpBeqv);  %AAB- must be zero
            d2Sbus_dshdpBeqv(kk,k) = (V.*conj(d2Ybus_dshdpBeqv*V)).'*lam;    %AAB- must be zero           
        end
        for kk=1:nPfsh                        
            %% Second Derivatives 
            d2Yff_dshdpsh = zeros(nl,1);                                       %AAB- must be zero
            d2Yft_dshdpsh = zeros(nl,1);                                       %AAB- must be zero
            d2Ytf_dshdpsh = zeros(nl,1);                                      %AAB- must be zero
            d2Ytt_dshdpsh = zeros(nl,1);                                     %AAB- must be zero
            
            d2Yf_dshdpsh = d2Yff_dshdpsh.*Cf + d2Yft_dshdpsh.*Ct;                
            d2Yt_dshdpsh = d2Ytf_dshdpsh.*Cf + d2Ytt_dshdpsh.*Ct;               
            
            d2Ybus_dshdpsh = Cf'*(d2Yf_dshdpsh) + Ct'*(d2Yt_dshdpsh);     
            d2Sbus_dshdpsh(kk,k) = (V.*conj(d2Ybus_dshdpsh*V)).'*lam;                  
        end
        for kk=1:nQtma
            QtmaPfdpSel=diagQtmaSel(:,iQtma(kk)).*diagPfdpSel(:,iPfdp(k)); %AAB- Selects only the active element controlling Droop and Qt with sh and ma
            YsmashdpSel=Ys.*QtmaPfdpSel;
            %% Second Derivatives 
            d2Yff_dshdpqtma = zeros(nl,1);                                       %AAB- must be zero
            d2Yft_dshdpqtma = ( 1j.*YsmashdpSel)./( k2.*abs(tap).*conj(tap) );   %AAB- Only active for elements with both Droop and Qt control
            d2Ytf_dshdpqtma = (-1j.*YsmashdpSel)./( k2.*abs(tap).*     tap  );  %AAB- Only active for elements with both Droop and Qt control
            d2Ytt_dshdpqtma = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dshdpqtma = d2Yff_dshdpqtma.*Cf + d2Yft_dshdpqtma.*Ct;                
            d2Yt_dshdpqtma = d2Ytf_dshdpqtma.*Cf + d2Ytt_dshdpqtma.*Ct;               
            
            d2Ybus_dvtmaqtma = Cf'*(d2Yf_dshdpqtma) + Ct'*(d2Yt_dshdpqtma);     
            d2Sbus_dshdpqtma(kk,k) = (V.*conj(d2Ybus_dvtmaqtma*V)).'*lam;                  
        end
        for kk=1:nVtma
            VtmaPfdpSel=diagVtmaSel(:,iVtma(kk)).*diagPfdpSel(:,iPfdp(k)); %AAB- Selects only the active element controlling Droop and Vt with sh and ma
            YsmashdpSel=Ys.*VtmaPfdpSel;
            
            %% Second Derivatives 
            d2Yff_dshdpvtma = zeros(nl,1);                                       %AAB- must be zero
            d2Yft_dshdpvtma = ( 1j.*YsmashdpSel)./( k2.*abs(tap).*conj(tap) );   %AAB- Only active for elements with both Droop and Qt control
            d2Ytf_dshdpvtma = (-1j.*YsmashdpSel)./( k2.*abs(tap).*     tap  );  %AAB- Only active for elements with both Droop and Qt control
            d2Ytt_dshdpvtma = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dshdpvtma = d2Yff_dshdpvtma.*Cf + d2Yft_dshdpvtma.*Ct;                
            d2Yt_dshdpvtma = d2Ytf_dshdpvtma.*Cf + d2Ytt_dshdpvtma.*Ct;               
            
            d2Ybus_dvtmavtma = Cf'*(d2Yf_dshdpvtma) + Ct'*(d2Yt_dshdpvtma);     
            d2Sbus_dshdpvtma(kk,k) = (V.*conj(d2Ybus_dvtmavtma*V)).'*lam;                  
        end
        for kk=1:nPfdp
            ShdpSel2=diagPfdpSel(:,iPfdp(kk)); %AAB- Selects the column of diagShSel representing only the active Sh
            Ysshdp=diagYsshdp(:,iPfdp(k)); %AAB- Selects the column of diagShsel representing only the active Sh
            
            %% Second Derivatives         
            d2Yff_dshdp2 = sparse( zeros(nl,1) );
            d2Yft_dshdp2 = sparse( (Ysshdp.*ShdpSel2)./(k2.*conj(tap)) );        %AAB- Only when Shifter selector 1 and Shifter selector 2 match there will be a derivative (When k = kk). Otherwise zero 
            d2Ytf_dshdp2 = sparse( (Ysshdp.*ShdpSel2)./(k2.*tap      ) );        %AAB- Only when Shifter selector 1 and Shifter selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Ytt_dshdp2 = sparse( zeros(nl,1) );
            
            d2Yf_dshdp2 = d2Yff_dshdp2.*Cf + d2Yft_dshdp2.*Ct;          
            d2Yt_dshdp2 = d2Ytf_dshdp2.*Cf + d2Ytt_dshdp2.*Ct;          
            
            d2Ybus_dshdp2 = Cf'*d2Yf_dshdp2 + Ct'*d2Yt_dshdp2;                   
            d2Sbus_dshdp2(kk,k) = (V.*conj(d2Ybus_dshdp2*V)).'*lam;            %AAB- size of [nVtma, nVtma]   
            
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
G110 = sparse(d2Sbus_dshdpVa);
G210 = sparse(d2Sbus_dshdpVm);
G510 = sparse(d2Sbus_dshdpBeqz);
G610 = sparse(d2Sbus_dshdpBeqv);
G710 = sparse(d2Sbus_dshdpsh);
G810 = sparse(d2Sbus_dshdpqtma);
G910 = sparse(d2Sbus_dshdpvtma);
G101 = G110.';
G102 = G210.';
G105 = G510.';
G106 = G610.';
G107 = G710.';
G108 = G810.';
G109 = G910.';
G1010 = sparse(d2Sbus_dshdp2);



