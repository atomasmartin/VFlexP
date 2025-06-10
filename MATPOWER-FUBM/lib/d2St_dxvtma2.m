function [H17, H27, H37, H47, H57, H67, H71, H72, H73, H74, H75, H76, H77] = d2St_dxvtma2(branch, V, mu, vcart)
%D2ST_DXVTMA2  Computes 2nd derivatives of complex brch power flow "from" w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma, vtmavtma.
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 4th argument. So far only Polar has been coded
%   
%   [HvtmaVa, HvtmaVm, HvtmaBz, HvtmaBv, HvtmaSh, Hvtmaqtma, HVavtma, HVmvtma, HBzvtma, HBvvtma, HShvtma, Hqtmavtma, Hvtmavtma] = D2ST_DXVTMA2(BRANCH, V, MU)
%   [HvtmaVa, HvtmaVm, HvtmaBz, HvtmaBv, HvtmaSh, Hvtmaqtma, HVavtma, HVmvtma, HBzvtma, HBvvtma, HShvtma, Hqtmavtma, Hvtmavtma] = D2ST_DXVTMA2(BRANCH, V, MU, 0)
%
%   Returns 13 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector MU with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HvtmaVa, HvtmaVm, HvtmaBz, HvtmaBv, HvtmaSh, Hvtmaqtma, HVavtma, HVmvtma, HBzvtma, HBvvtma, HShvtma, Hqtmavtma, Hvtmavtma] = D2ST_DXVTMA2(BRANCH, V, MU, 1)
%
%   Not Coded Yet
%
%   Examples:
%       il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
%       [HvtmaVa, HvtmaVm, HvtmaBz, HvtmaBv, HvtmaSh, Hvtmaqtma, HVavtma, ...
%                 HVmvtma, HBzvtma, HBvvtma, HShvtma, Hqtmavtma, Hvtmavtma] = ...
%                                                                    d2St_dxvtma2(branch(il,:), V, mu, 0);
%
%       Here the output matrices correspond to:
%           HvtmaVa   = d/dvtma (dSt_dVa.'   * mu)
%           HvtmaVm   = d/dvtma (dSt_dVm.'   * mu)
%           HvtmaBz   = d/dvtma (dSt_dBeqz.' * mu)
%           HvtmaBv   = d/dvtma (dSt_dBeqv.' * mu)
%           HvtmaSh   = d/dvtma (dSt_dSh.'   * mu)
%           Hvtmaqtma = d/dvtma (dSt_dqtma.' * mu)
%           HVavtma   = d/dVa   (dSt_dvtma.' * mu)
%           HVmvtma   = d/dVm   (dSt_dvtma.' * mu)
%           HBzvtma   = d/dBeqz (dSt_dvtma.' * mu)
%           HBvvtma   = d/dBeqv (dSt_dvtma.' * mu)
%           HShvtma   = d/dSh   (dSt_dvtma.' * mu)
%           Hqtmavtma = d/dqtma (dSt_dvtma.' * mu)
%           Hvtmavtma = d/dvtma (dSt_dvtma.' * mu)
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
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2St_dxvtma2: Derivatives of Flow Limit equations w.r.t ma in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
     
    %Selector of active ma/tap 
    QtmaSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector ma controlling Qt
    VtmaSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector ma controlling Vt
    BeqzSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqz controlling Zero Constraint
    BeqvSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqv controlling Vf
    PfshSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector sh controlling Pf  
    
    QtmaSel(iQtma) = 1;                        %AAB- Fill the selector with 1 where ma is active and controlling Qt
    VtmaSel(iVtma) = 1;                        %AAB- Fill the selector with 1 where ma is active and controlling Vt
    BeqzSel(iBeqz) = 1;                        %AAB- Fill the selector with 1 where Beqz is active and controlling Zero Constraint
    BeqvSel(iBeqv) = 1;                        %AAB- Fill the selector with 1 where Beqv is active and controlling Vf
    PfshSel(iPfsh) = 1;                        %AAB- Fill the selector with 1 where sh is active and controlling Pf
    
    diagQtmaSel = sparse( diag(QtmaSel));      %AAB- Diagonal of the selector where ma is active and controlling Qt
    diagVtmaSel = sparse( diag(VtmaSel));      %AAB- Diagonal of the selector where ma is active and controlling Qt
    diagPfshSel = sparse( diag(PfshSel));      %AAB- Diagonal of the selector where sh is active and controlling Pf
    diagBeqzSel = sparse( diag(BeqzSel));      %AAB- Diagonal of the selector where Beqz is active and controlling Zero Constraint
    diagBeqvSel = sparse( diag(BeqvSel));      %AAB- Diagonal of the selector where Beqv is active and controlling Vf   
    
    diagYsVtma   = sparse( diag(VtmaSel.*Ys) );        %AAB- ma/tap selector multilied by the series addmitance Ys, size [nl,nl]
    diagYttBVtma = sparse( diag(VtmaSel.*YttB) );      %AAB- ma selector multilied by the series addmitance Ytt, size [nl,nl]

    %Dimensionalize (Allocate for computational speed)
    dYff_dvtma = sparse( zeros(nl,nVtma) );
    dYft_dvtma = sparse( zeros(nl,nVtma) );
    dYtf_dvtma = sparse( zeros(nl,nVtma) ); 
    dYtt_dvtma = sparse( zeros(nl,nVtma) );
    
    d2St_dvtmaVa   = sparse( zeros(nb,   nVtma) ); 
    d2St_dvtmaVm   = sparse( zeros(nb,   nVtma) );
    d2St_dvtmaBeqz = sparse( zeros(nBeqz,nVtma) );
    d2St_dvtmaBeqv = sparse( zeros(nBeqv,nVtma) );
    d2St_dvtmash   = sparse( zeros(nPfsh,nVtma) );   
    d2St_dvtmaqtma = sparse( zeros(nQtma,nVtma) );  
    d2St_dvtma2    = sparse( zeros(nVtma,nVtma) );
    
    for k=1:nVtma 
        for kk=1:nb %dQtmaVx
            %% Second Derivatives
            Ysvtma=diagYsVtma(:,iVtma(k)); %AAB- Selects the column of diagShsel representing only the active Sh
            YttBvtma=diagYttBVtma(:,iVtma(k)); %AAB- Selects the column of diagmaAux representing only the active ma for the specified control
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. ma
            %dYff_dvtma(:, k) = sparse( -2*YttBvtma./( (k2.^2).*((abs(tap)).^3) ) );
            %dYft_dvtma(:, k) = sparse( Ysvtma./( k2.*(abs(tap).*conj(tap)) ) ); 
            dYtf_dvtma(:, k) = sparse( Ysvtma./( k2.*(abs(tap).*     tap ) ) ); 
            dYtt_dvtma(:, k) = sparse( zeros(nl,1) );
        
            %Partials of Yf, Yt w.r.t. ma
            %dYf_dvtma = dYff_dvtma(:, k).* Cf + dYft_dvtma(:, k).* Ct; %AAB- size [nl,nb] per active ma
            dYt_dvtma = dYtf_dvtma(:, k).* Cf + dYtt_dvtma(:, k).* Ct; %AAB- size [nl,nb] per active ma

            %2nd Derivatives of Sf w.r.t. VtmaVx
            d2St_dvtmaVa(kk, k) = ((diag(Ct*V)*conj(dYt_dvtma*dVa(:,kk))   +   (Ct*dVa(:,kk)).*conj(dYt_dvtma*V)).')*mu; %AAB- Final dSt_dmaVa has a size of [nl, nVtma] 
            d2St_dvtmaVm(kk, k) = ((diag(Ct*V)*conj(dYt_dvtma*dVm(:,kk))   +   (Ct*dVm(:,kk)).*conj(dYt_dvtma*V)).')*mu; %AAB- Final dSt_dmaVm has a size of [nl, nVtma]  
        end      
        for kk=1:nBeqz
            %VtmaBeqzSel=diagVtmaSel(:,iVtma(k)).*diagBeqzSel(:,iBeqz(kk)); %AAB- Selects the column of diagVtmaBeqzSel representing only the active element controlling Zero Constraint and Vt with Beqz and ma
            
            %% Second Derivatives %The vtmaBeqz derivative is zero because Yft, Ytf and Ytt do not share ma and Beqz for any case.
            
            %d2Yff_dvtmaBeqz = (-2j*VtmaBeqzSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dvtmaBeqz = zeros(nl,1);
            %d2Yft_dvtmaBeqz = zeros(nl,1);                                  %AAB- must be zero
            d2Ytf_dvtmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dvtmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            
            %d2Yf_dvtmaBeqz = d2Yff_dvtmaBeqz.*Cf + d2Yft_dvtmaBeqz.*Ct;     %AAB- must be zero
            d2Yt_dvtmaBeqz = d2Ytf_dvtmaBeqz.*Cf + d2Ytt_dvtmaBeqz.*Ct;      %AAB- must be zero
            
            d2St_dvtmaBeqz(kk,k) = (diag(Ct*V)*conj(d2Yt_dvtmaBeqz*V)).'*mu; %AAB- must be zero
        end
        for kk=1:nBeqv
            %VtmaBeqvSel=diagVtmaSel(:,iVtma(k)).*diagBeqvSel(:,iBeqv(kk)); %AAB- Selects the column of diagVtmaBeqvSel representing only the active element controlling Vf and Vt with Beqv and ma
           
            %% Second Derivatives %The vtmaBeqv derivative is zero becauseYff, Yft, Ytf and Ytt do not share ma and Beqv for any case.
            %d2Yff_dvtmaBeqv = (-2j*VtmaBeqvSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dvtmaBeqv = zeros(nl,1);
            %d2Yft_dvtmaBeqv = zeros(nl,1);                                  %AAB- must be zero
            d2Ytf_dvtmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dvtmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            
            %d2Yf_dvtmaBeqv = d2Yff_dvtmaBeqv.*Cf + d2Yft_dvtmaBeqv.*Ct;     %AAB- must be zero
            d2Yt_dvtmaBeqv = d2Ytf_dvtmaBeqv.*Cf + d2Ytt_dvtmaBeqv.*Ct;      %AAB- must be zero
            
            d2St_dvtmaBeqv(kk,k) = (diag(Ct*V)*conj(d2Yt_dvtmaBeqv*V)).'*mu; %AAB- must be zero
        end
        for kk=1:nPfsh
            VtmaPfshSel=diagVtmaSel(:,iVtma(k)).*diagPfshSel(:,iPfsh(kk)); %AAB- Selects only the active element controlling Pf and Vt with sh and ma
            YsVtmashSel=Ys.*VtmaPfshSel;
            
            %% Second Derivatives 
            %d2Yff_dvtmash = zeros(nl,1);                                       %AAB- must be zero
            %d2Yft_dvtmash = ( 1j.*YsvtmashSel)./( k2.*abs(tap).*conj(tap) );  %AAB- Only active for elements with both Pf and Qt control
            d2Ytf_dvtmash = (-1j.*YsVtmashSel)./( k2.*abs(tap).*     tap  );   %AAB- Only active for elements with both Pf and Qt control
            d2Ytt_dvtmash = zeros(nl,1);                                        %AAB- must be zero
            
            %d2Yf_dvtmash = d2Yff_dvtmash.*Cf + d2Yft_dvtmash.*Ct;                
            d2Yt_dvtmash = d2Ytf_dvtmash.*Cf + d2Ytt_dvtmash.*Ct;               
            
            d2St_dvtmash(kk,k) = (diag(Ct*V)*conj(d2Yt_dvtmash*V)).'*mu;                  
        end
        for kk=1:nQtma
            %% Second Derivatives 
            %d2Yff_dvtmaqtma = zeros(nl,1);                                %AAB- must be zero
            %d2Yft_dvtmaqtma = zeros(nl,1);                                %AAB- must be zero
            d2Ytf_dvtmaqtma = zeros(nl,1);                                 %AAB- must be zero
            d2Ytt_dvtmaqtma = zeros(nl,1);                                 %AAB- must be zero
            
            %d2Yf_dvtmaqtma = d2Yff_dvtmaqtma.*Cf + d2Yft_dvtmaqtma.*Ct;                
            d2Yt_dvtmaqtma = d2Ytf_dvtmaqtma.*Cf + d2Ytt_dvtmaqtma.*Ct;               
            
            d2St_dvtmaqtma(kk,k) = (diag(Ct*V)*conj(d2Yt_dvtmaqtma*V)).'*mu;                  
        end
        for kk=1:nVtma
            vtmaSel2=diagVtmaSel(:,iVtma(kk)); %AAB- Selects the column of diagmaSel representing only the active ma
            
            %% Second Derivatives         
            %d2Yff_dma2 = sparse( ( 6.*YttBvtma.*vtmaSel2 )./( (k2.^2).*((abs(tap)).^4) ) );       %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            %d2Yft_dma2 = sparse( ( -2.*Ysvtma.*vtmaSel2 )./( k2.*((abs(tap)).^2).*conj(tap) ) );  %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero 
            d2Ytf_dma2 = sparse( ( -2.*Ysvtma.*vtmaSel2 )./( k2.*((abs(tap)).^2).*     tap  ) ); %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Ytt_dma2 = sparse( zeros(nl,1) );
            
            %d2Yf_dma2 = d2Yff_dma2.*Cf + d2Yft_dma2.*Ct;          
            d2Yt_dma2 = d2Ytf_dma2.*Cf + d2Ytt_dma2.*Ct;          
                  
            d2St_dvtma2(kk,k) = (diag(Ct*V)*conj(d2Yt_dma2*V)).'*mu;         %AAB- size of [nVtma, nVtma]   
            
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
H17 = sparse(d2St_dvtmaVa);
H27 = sparse(d2St_dvtmaVm);
H37 = sparse(d2St_dvtmaBeqz);
H47 = sparse(d2St_dvtmaBeqv);
H57 = sparse(d2St_dvtmash);
H67 = sparse(d2St_dvtmaqtma);
H71 = H17.';
H72 = H27.';
H73 = H37.';
H74 = H47.';
H75 = H57.';
H76 = H67.';
H77 = sparse(d2St_dvtma2);











