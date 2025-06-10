function [G18, G28, G58, G68, G78, G81, G82, G85, G86, G87, G88] = d2Sbus_dxqtma2(branch, V, lam, vcart)
%D2SBUS_DXQTMA2   Computes 2nd derivatives of power injection w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma
%
%   The derivatives will be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument. So far only polar
%
%   [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = D2SBUS_DXQTMA2(BRANCH, V, LAM)

%   [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = D2SBUS_DXQTMA2(BRANCH, V, LAM, 0)
%
%   Returns 11 matrices containing the partial derivatives the product of a
%   vector LAM with the 1st partial derivatives of the complex  power 
%   injections w.r.t. GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, 
%   GBeqzqtma, GBeqvqtma, GShqtma and Gqtmaqtma
%
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [GqtmaVa, GqtmaVm, GqtmaBeqz, GqtmaBeqv, GqtmaSh, GVaqtma, GVmqtma, GBeqzqtma, GBeqvqtma, GShqtma, Gqtmaqtma] = ...
%                                                                           d2Sbus_dxqtma2(branch, V, lam, vcart);
%
%       Here the output matrices correspond to:
%           GqtmaVa = d/dVa   (dSbus_dqtma.' * lam)
%           GqtmaVm = d/dVm   (dSbus_dqtma.' * lam)
%           GqtmaBz = d/dBeqz (dSbus_dqtma.' * lam)
%           GqtmaBv = d/dBeqv (dSbus_dqtma.' * lam)
%           GqtmaSh = d/dsh   (dSbus_dqtma.' * lam)
%           GVaqtma = d/dqtma (dSbus_dVa.'   * lam)
%           GVmqtma = d/dqtma (dSbus_dVm.'   * lam)
%           GBzqtma = d/dqtma (dSbus_dBeqz.' * lam)
%           GBvqtma = d/dqtma (dSbus_dBeqv.' * lam)
%           GShqtma = d/dqtma (dSbus_dsh.'   * lam)
%           Gmaqtma = d/dqtma (dSbus_dqtma.' * lam)
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

[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% Calculation of derivatives
if vcart
    error('d2Sbus_dxqtma2: Derivatives of Power balance equations w.r.t ma/tap in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
    
    %Selector of active ma/tap 
    QtmaSel = sparse( zeros(nl,1) );             %AAB- Vector of zeros for the selector ma
    PfshSel = sparse( zeros(nl,1) );             %AAB- Vector of zeros for the selector sh    
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
    
    d2Sbus_dqtmaVa   = sparse( zeros(nb,   nQtma) ); 
    d2Sbus_dqtmaVm   = sparse( zeros(nb,   nQtma) ); 
    d2Sbus_dqtmaBeqz = sparse( zeros(nBeqz,nQtma) );
    d2Sbus_dqtmaBeqv = sparse( zeros(nBeqv,nQtma) );
    d2Sbus_dqtmash   = sparse( zeros(nPfsh,nQtma) );
    d2Sbus_dqtma2    = sparse( zeros(nQtma,nQtma) );
    
    for k=1:nQtma 
        for kk=1:nb %dQtmaVx
            %% Second Derivatives  
            Ysma=diagYsQtma(:,iQtma(k)); %AAB- Selects the column of diagmasel representing only the active ma
            YttBma=diagYttBQtma(:,iQtma(k)); %AAB- Selects the column of diagmaAux representing only the active ma for the specified control
        
            %Partials of Ytt, Yff, Yft and Ytf w.r.t. ma
            dYff_dma(:, k) = sparse( -2*YttBma./( (k2.^2).*((abs(tap)).^3) ) );
            dYft_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*conj(tap)) ) ); 
            dYtf_dma(:, k) = sparse( Ysma./( k2.*(abs(tap).*     tap ) ) ); 
            dYtt_dma(:, k) = sparse( zeros(nl,1) );
        
            %Partials of Yf, Yt, Ybus w.r.t. ma
            dYf_dma = dYff_dma(:, k).* Cf + dYft_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma
            dYt_dma = dYtf_dma(:, k).* Cf + dYtt_dma(:, k).* Ct; %AAB- size [nl,nb] per active ma

            dYbus_dma = Cf' * dYf_dma + Ct' * dYt_dma;     %AAB- size [nb,nb] per active ma       

            %2nd Derivatives of Sbus w.r.t. qtmaVx
            d2Sbus_dqtmaVa(kk, k) = ((dVa(:,kk).*conj(dYbus_dma*V) + V.*conj(dYbus_dma*(dVa(:,kk)))).')*lam;       %AAB- Final d2Sbus_dmaVa has a size of [nb, nQtma] 
            d2Sbus_dqtmaVm(kk, k) = ((dVm(:,kk).*conj(dYbus_dma*V) + V.*conj(dYbus_dma*(dVm(:,kk)))).')*lam; %AAB- Final d2Sbus_dmaVm has a size of [nb, nQtma] 
        end
        for kk=1:nBeqz
            qtmaBeqzSel=diagQtmaSel(:,iQtma(k)).*diagBeqzSel(:,iBeqz(kk)); %AAB- Selects only the active element controlling Zero Constraint and Qt with Beq and ma

            %% Second Derivatives %The qtmaBeqz derivative is zero because Yft, Ytf and Ytt do not share ma and Beqz for any case.
            d2Yff_dmaBeqz = (-2j*qtmaBeqzSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dmaBeqz = zeros(nl,1);
            d2Yft_dmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytf_dmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dmaBeqz = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dmaBeqz = d2Yff_dmaBeqz.*Cf + d2Yft_dmaBeqz.*Ct;          
            d2Yt_dmaBeqz = d2Ytf_dmaBeqz.*Cf + d2Ytt_dmaBeqz.*Ct;          %AAB- must be zero
            
            d2Ybus_dmaBeqz = Cf'*(d2Yf_dmaBeqz) + Ct'*(d2Yt_dmaBeqz);      
            d2Sbus_dqtmaBeqz(kk,k) = (V.*conj(d2Ybus_dmaBeqz*V)).'*lam;              
        end
        for kk=1:nBeqv
            qtmaBeqvSel=diagQtmaSel(:,iQtma(k)).*diagBeqvSel(:,iBeqv(kk)); %AAB- Selects only the active element controlling Vf and Qt with Beq and ma

            %% Second Derivatives %The qtmaBeqv derivative is zero because Yft, Ytf and Ytt do not share ma and Beqv for any case.
            d2Yff_dmaBeqv = (-2j*qtmaBeqvSel)./( (k2.^2).*((abs(tap)).^3) );    %AAB-%Original: d2Yff_dmaBeqv = zeros(nl,1);
            d2Yft_dmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytf_dmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            d2Ytt_dmaBeqv = zeros(nl,1);                                   %AAB- must be zero
            
            d2Yf_dmaBeqv = d2Yff_dmaBeqv.*Cf + d2Yft_dmaBeqv.*Ct;          
            d2Yt_dmaBeqv = d2Ytf_dmaBeqv.*Cf + d2Ytt_dmaBeqv.*Ct;          %AAB- must be zero
            
            d2Ybus_dmaBeqv = Cf'*(d2Yf_dmaBeqv) + Ct'*(d2Yt_dmaBeqv);      
            d2Sbus_dqtmaBeqv(kk,k) = (V.*conj(d2Ybus_dmaBeqv*V)).'*lam;               
        end
        for kk=1:nPfsh
            qtmaPfshSel=diagQtmaSel(:,iQtma(k)).*diagPfshSel(:,iPfsh(kk)); %AAB- Selects only the active element controlling Pf and Qt with sh and ma
            YsmashSel2=Ys.*qtmaPfshSel;                                    %AAB- Selector multiplied by Ys
            
            %% Second Derivatives 
            d2Yff_dmash = zeros(nl,1);                                     %AAB- must be zero
            d2Yft_dmash = ( 1j.*YsmashSel2)./( k2.*abs(tap).*conj(tap) );  %AAB- Only active for elements with both Pf and Qt control
            d2Ytf_dmash = (-1j.*YsmashSel2)./( k2.*abs(tap).*     tap  );  %AAB- Only active for elements with both Pf and Qt control
            d2Ytt_dmash = zeros(nl,1);                                     %AAB- must be zero
            
            d2Yf_dmash = d2Yff_dmash.*Cf + d2Yft_dmash.*Ct;                
            d2Yt_dmash = d2Ytf_dmash.*Cf + d2Ytt_dmash.*Ct;               
            
            d2Ybus_dmash = Cf'*(d2Yf_dmash) + Ct'*(d2Yt_dmash);     
            d2Sbus_dqtmash(kk,k) = (V.*conj(d2Ybus_dmash*V)).'*lam;                  
        end
        for kk=1:nQtma
            maSel2=diagQtmaSel(:,iQtma(kk)); %AAB- Selects the column of diagmaSel representing only the active ma
            
            %% Second Derivatives         
            d2Yff_dma2 = sparse( ( 6.*YttBma.*maSel2 )./( (k2.^2).*((abs(tap)).^4) ) );      %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Yft_dma2 = sparse( ( -2.*Ysma.*maSel2 )./( k2.*((abs(tap)).^2).*conj(tap) ) ); %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero 
            d2Ytf_dma2 = sparse( ( -2.*Ysma.*maSel2 )./( k2.*((abs(tap)).^2).*     tap  ) ); %AAB- Only when ma selector 1 and ma selector 2 match there will be a derivative (When k = kk). Otherwise zero
            d2Ytt_dma2 = sparse( zeros(nl,1) );
            
            d2Yf_dma2 = d2Yff_dma2.*Cf + d2Yft_dma2.*Ct;          
            d2Yt_dma2 = d2Ytf_dma2.*Cf + d2Ytt_dma2.*Ct;          
            
            d2Ybus_dma2 = Cf'*d2Yf_dma2 + Ct'*d2Yt_dma2;                   
            d2Sbus_dqtma2(kk,k) = (V.*conj(d2Ybus_dma2*V)).'*lam;            %AAB- size of [nQtma, nQtma]   
            
        end
    end
end

%Assigning the partial derivatives with their respective outputs. 
G18 = sparse(d2Sbus_dqtmaVa);
G28 = sparse(d2Sbus_dqtmaVm);
G58 = sparse(d2Sbus_dqtmaBeqz);
G68 = sparse(d2Sbus_dqtmaBeqv);
G78 = sparse(d2Sbus_dqtmash);
G81 = G18.';
G82 = G28.';
G85 = G58.';
G86 = G68.';
G87 = G78.';
G88 = sparse(d2Sbus_dqtma2);