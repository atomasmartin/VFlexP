function [H16, H26, H36, H46, H56, H61, H62, H63, H64, H65, H66] = d2Sf_dxqtma2(branch, V, mu, vcart)
%D2SF_DXQTMA2  Computes 2nd derivatives of complex brch power flow "from" w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma.
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 4th argument. So far only Polar has been coded
%   
%   [HqtmaVa, HqtmaVm, HqtmaBz, HqtmaBv, HqtmaSh, HVaqtma, HVmqtma, HBzqtma, HBvqtma, HShqtma, Hqtmaqtma] = D2SF_DXQTMA2(BRANCH, V, MU)
%   [HqtmaVa, HqtmaVm, HqtmaBz, HqtmaBv, HqtmaSh, HVaqtma, HVmqtma, HBzqtma, HBvqtma, HShqtma, Hqtmaqtma] = D2SF_DXQTMA2(BRANCH, V, MU, 0)
%
%   Returns 11 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector MU with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HqtmaVa, HqtmaVm, HqtmaBz, HqtmaBv, HqtmaSh, HVaqtma, HVmqtma, HBzqtma, HBvqtma, HShqtma, Hqtmaqtma] = D2SF_DXQTMA2(BRANCH, V, MU, 1)
%
%   Not Coded Yet
%
%   Examples:
%       il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
%       [HqtmaVa, HqtmaVm, HqtmaBz, HqtmaBv, HqtmaSh, HVaqtma, HVmqtma, HBzqtma, HBvqtma, HShqtma, Hqtmaqtma] = d2Sf_dxqtma2(branch(il,:), V, mu, 0);
%
%       Here the output matrices correspond to:
%           HqtmaVa   = d/dqtma (dSf_dVa.'   * mu)
%           HqtmaVm   = d/dqtma (dSf_dVm.'   * mu)
%           HqtmaBz   = d/dqtma (dSf_dBeqz.' * mu)
%           HqtmaBv   = d/dqtma (dSf_dBeqv.' * mu)
%           HqtmaSh   = d/dqtma (dSf_dSh.'   * mu)
%           HVaqtma   = d/dVa   (dSf_dqtma.' * mu)
%           HVmqtma   = d/dVm   (dSf_dqtma.' * mu)
%           HBzqtma   = d/dBeqz (dSf_dqtma.' * mu)
%           HBvqtma   = d/dBeqv (dSf_dqtma.' * mu)
%           HShqtma   = d/dSh   (dSf_dqtma.' * mu)
%           Hqtmaqtma = d/dqtma (dSf_dqtma.' * mu)
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
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sf_dxma2: Derivatives of Flow Limit equations w.r.t ma in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active ma/tap 
    QtmaSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector ma
    PfshSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector sh    
    BeqzSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqz  
    BeqvSel = sparse( zeros(nl,1) );           %AAB- Vector of zeros for the selector Beqv    
    
    QtmaSel(iQtma) = 1;                        %AAB- Fill the selector with 1 where ma is active and controlling Qt
    PfshSel(iPfsh) = 1;                        %AAB- Fill the selector with 1 where sh is active and controlling Pf
    BeqzSel(iBeqz) = 1;                        %AAB- Fill the selector with 1 where Beqz is active and controlling Zero Constraint
    BeqvSel(iBeqv) = 1;                        %AAB- Fill the selector with 1 where Beqv is active and controlling Vf
    
    diagQtmaSel = sparse( diag(QtmaSel));      %AAB- Diagonal of the selector where ma is active and controlling Qt
    diagPfshSel = sparse( diag(PfshSel));      %AAB- Diagonal of the selector where sh is active and controlling Pf
    diagBeqzSel = sparse( diag(BeqzSel));      %AAB- Diagonal of the selector where Beqz is active and controlling Zero Constraint
    diagBeqvSel = sparse( diag(BeqvSel));      %AAB- Diagonal of the selector where Beqv is active and controlling Vf   
    
    diagYsQtma = sparse( diag(QtmaSel.*Ys) );      %AAB- ma/tap selector multilied by the series addmitance Ys, size [nl,nl]
    diagYttBQtma= sparse( diag(QtmaSel.*YttB) );   %AAB- ma selector multilied by the series addmitance Ytt, size [nl,nl]
    
    %Dimensionalize (Allocate for computational speed)
    dYtt_dma = sparse( zeros(nl,nQtma) );
    dYff_dma = sparse( zeros(nl,nQtma) );
    dYft_dma = sparse( zeros(nl,nQtma) );
    dYtf_dma = sparse( zeros(nl,nQtma) ); 
    
    d2Sf_dqtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Sf_dqtmaVm   = sparse( zeros(nb,   nQtma) );
    d2Sf_dqtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Sf_dqtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Sf_dqtmash   = sparse( zeros(nPfsh,nQtma) );    
    d2Sf_dqtma2    = sparse( zeros(nQtma,nQtma) );
    
    for k=1:nQtma 
        for kk=1:nb %dQtmaVx
            %% Second Derivatives 
            Ysma=diagYsQtma(:,iQtma(k)); %AAB- Selects the column of diagShsel representing only the active Sh
            YttBma=diagYttBQtma(:,iQtma(k)); %AAB- Selects the column of diagmaAux representing only the active ma for the specified control
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. ma
            dYff_dma(:, k) = sparse( -2*YttBma./( (k2.^2).*((abs(tap)).^3) ) );
            dYft_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*conj(tap)) ) ); 
            %dYtf_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*     tap ) ) ); 
            %dYtt_dma(:, k) = sparse( zeros(nl,1) );
        
            %Partials of Yf, Yt w.r.t. ma
            dYf_dma = dYff_dma(:, k).* Cf + dYft_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma
            %dYt_dma = dYtf_dma(:, k).* Cf + dYtt_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma

            %2nd Derivatives of Sf w.r.t. qtmaVx
            d2Sf_dqtmaVa(kk, k) = ((diag(Cf*V)*conj(dYf_dma*dVa(:,kk))   +   (Cf*dVa(:,kk)).*conj(dYf_dma*V)).')*mu; %AAB- Final dSf_dmaVa has a size of [nl, nQtma] 
            d2Sf_dqtmaVm(kk, k) = ((diag(Cf*V)*conj(dYf_dma*dVm(:,kk))   +   (Cf*dVm(:,kk)).*conj(dYf_dma*V)).')*mu; %AAB- Final dSf_dmaVm has a size of [nl, nQtma] 
        end
        for kk=1:nBeqz
            QtmaBeqzSel=diagQtmaSel(:,iQtma(k)).*diagBeqzSel(:,iBeqz(kk)); %AAB- Selects the column of diagqtmaBeqzSel representing only the active element controlling Zero Constraint and Qt with Beqz and ma
            
            %% Second Derivatives %The qtmaBeqz derivative is zero because Yft, Ytf and Ytt do not share ma and Beqz for any case.
            d2Yff_dmaBeqz = (-2j*QtmaBeqzSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dmaBeqz = zeros(nl,1);
            d2Yft_dmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            %d2Ytf_dmaBeqz = zeros(nl,1);                                  %AAB- must be zero
            %d2Ytt_dmaBeqz = zeros(nl,1);                                  %AAB- must be zero
            
            d2Yf_dmaBeqz = d2Yff_dmaBeqz.*Cf + d2Yft_dmaBeqz.*Ct;          %AAB- must be zero
            %d2Yt_dmaBeqz = d2Ytf_dmaBeqz.*Cf + d2Ytt_dmaBeqz.*Ct;         %AAB- must be zero
            
            d2Sf_dqtmaBeqz(kk,k) = (diag(Cf*V)*conj(d2Yf_dmaBeqz*V)).'*mu;   %AAB- must be zero
        end
        for kk=1:nBeqv
            QtmaBeqvSel=diagQtmaSel(:,iQtma(k)).*diagBeqvSel(:,iBeqv(kk)); %AAB- Selects the column of diagqtmaBeqvSel representing only the active element controlling Vf and Qt with Beqv and ma
            
            %% Second Derivatives %The qtmaBeqv derivative is zero because Yft, Ytf and Ytt do not share ma and Beqv for any case.
            d2Yff_dmaBeqv = (-2j*QtmaBeqvSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dmaBeqv = zeros(nl,1);
            d2Yft_dmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            %d2Ytf_dmaBeqv = zeros(nl,1);                                  %AAB- must be zero
            %d2Ytt_dmaBeqv = zeros(nl,1);                                  %AAB- must be zero
            
            d2Yf_dmaBeqv = d2Yff_dmaBeqv.*Cf + d2Yft_dmaBeqv.*Ct;          %AAB- must be zero
            %d2Yt_dmaBeqv = d2Ytf_dmaBeqv.*Cf + d2Ytt_dmaBeqv.*Ct;         %AAB- must be zero
            
            d2Sf_dqtmaBeqv(kk,k) = (diag(Cf*V)*conj(d2Yf_dmaBeqv*V)).'*mu;   %AAB- must be zero
        end
        for kk=1:nPfsh
            QtmaPfshSel=diagQtmaSel(:,iQtma(k)).*diagPfshSel(:,iPfsh(kk)); %AAB- Selects only the active element controlling Pf and Qt with sh and ma
            YsmashSel=Ys.*QtmaPfshSel;
            
            %% Second Derivatives 
            d2Yff_dmash = zeros(nl,1);                                     %AAB- must be zero
            d2Yft_dmash = ( 1j.*YsmashSel)./( k2.*abs(tap).*conj(tap) );   %AAB- Only active for elements with both Pf and Qt control
            %d2Ytf_dmash = (-1j.*YsmashSel)./( k2.*abs(tap).*     tap  );  %AAB- Only active for elements with both Pf and Qt control
            %d2Ytt_dmash = zeros(nl,1);                                    %AAB- must be zero
            
            d2Yf_dmash = d2Yff_dmash.*Cf + d2Yft_dmash.*Ct;                
            %d2Yt_dmash = d2Ytf_dmash.*Cf + d2Ytt_dmash.*Ct;               
            
            d2Sf_dqtmash(kk,k) = (diag(Cf*V)*conj(d2Yf_dmash*V)).'*mu;                  
        end
        for kk=1:nQtma
            maSel2=diagQtmaSel(:,iQtma(kk)); %AAB- Selects the column of diagmaSel representing only the active ma
            
            %% Second Derivatives         
            d2Yff_dma2 = sparse( ( 6.*YttBma.*maSel2 )./( (k2.^2).*((abs(tap)).^4) ) );       %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Yft_dma2 = sparse( ( -2.*Ysma.*maSel2 )./( k2.*((abs(tap)).^2).*conj(tap) ) );  %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero 
            %d2Ytf_dma2 = sparse( ( -2.*Ysma.*maSel2 )./( k2.*((abs(tap)).^2).*     tap  ) ); %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            %d2Ytt_dma2 = sparse( zeros(nl,1) );
            
            d2Yf_dma2 = d2Yff_dma2.*Cf + d2Yft_dma2.*Ct;          
            %d2Yt_dma2 = d2Ytf_dma2.*Cf + d2Ytt_dma2.*Ct;          
                  
            d2Sf_dqtma2(kk,k) = (diag(Cf*V)*conj(d2Yf_dma2*V)).'*mu;         %AAB- size of [nQtma, nQtma]   
            
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
H16 = sparse(d2Sf_dqtmaVa);
H26 = sparse(d2Sf_dqtmaVm);
H36 = sparse(d2Sf_dqtmaBeqz);
H46 = sparse(d2Sf_dqtmaBeqv);
H56 = sparse(d2Sf_dqtmash);
H61 = H16.';
H62 = H26.';
H63 = H36.';
H64 = H46.';
H65 = H56.';
H66 = sparse(d2Sf_dqtma2);











