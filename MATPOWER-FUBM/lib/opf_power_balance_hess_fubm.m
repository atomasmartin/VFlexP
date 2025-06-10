function d2G = opf_power_balance_hess_fubm(x, lambda, mpc, mpopt)
%OPF_POWER_BALANCE_HESS_FUBM  Evaluates Hessian of power balance constraints.
%   D2G = OPF_POWER_BALANCE_HESS_AAB(X, LAMBDA, OM, MPOPT)
%
%   Hessian evaluation function for AC/DC active and reactive power balance
%   constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints (Pmis, Qmis, VfBeqmis)
%     MPC : MATPOWER case struct
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of power balance constraints.
%
%   Example:
%       d2G = opf_power_balance_hess(x, lambda, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_FCN, OPF_POWER_BALANCE_HESS.
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM-extra fields for FUBM
%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% Identifier of AC/DC grids
%%AAB--------------------------------------------------------------------- 
%%identifier of Zero Constraint VSCs
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with fixed Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
%identifier of converters with Ploss correction
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses

%%------------------------------------------------------------------------- 
%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift
%%------------------------------------------------------------------------- 

%% Reconstruction of V
if mpopt.opf.v_cartesian
    error('opf_power_balance_hess_aab: FUBM formulation with voltages in cartesian coordinates has not been coded yet.')
    %[Vr, Vi, Pg, Qg] = deal(x{:}); %AAB-This is not ready for FUBM
    %V = Vr + 1j * Vi;           %% reconstruct V
else %AAB- Polar variables
    %%AAB------------------------------------------------------------------
    [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, 1); %AAB- Deals optimisation variables
    V = Vm .* exp(1j * Va);     %% reconstruct V
    %%---------------------------------------------------------------------
end
%% AAB---------------------------------------------------------------------
%%update mpc.branch with FUBM from x
if nBeqz % AC/DC Formulation
    branch(iBeqz,BEQ)=Beqz; %AAB- Update the data from Beqz to the branch matrix
end
if nBeqv
    branch(iBeqv,BEQ)=Beqv; %AAB- Update the data from Beqv to the branch matrix  
end
if nPfsh
    branch(iPfsh,SHIFT) = ShAng*180/pi;  %AAB- Update the data from Theta_shift to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
end
if nQtma
    branch(iQtma,TAP) = maQt;  %AAB- Update the data from ma/tap to the branch matrix.
end
if nVtma
    branch(iVtma,TAP) = maVt;  %AAB- Update the data from ma/tap to the branch matrix.
end
if nPfdp
    branch(iPfdp,SHIFT) = ShAngDp*180/pi;  %AAB- Update the data from Theta_shift Droop to the branch matrix (It is returnded to degrees since inside makeYbus_aab it is converted to radians).
end
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %<<AAB-Ybus calculation with updated variables- Original: makeYbus
%% Standard IEC 62751-2 Ploss Correction for VSC losses
if nVscL
    %%compute branch power flows
    brf=branch(:, F_BUS);              %%AAB- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
    It= Yt(:, :) * V;                  %%AAB- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
    %%compute VSC Power Loss
    PLoss_IEC = branch(iVscL,ALPH3).*((abs(It(iVscL))).^2) + branch(iVscL,ALPH2).*((abs(It(iVscL)))) + branch(iVscL,ALPH1); %%AAB- Standard IEC 62751-2 Ploss Correction for VSC losses 
    branch(iVscL,GSW) = PLoss_IEC./(abs(V(brf(iVscL))).^2);    %%AAB- VSC Gsw Update
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch); %<<AAB-Ybus calculation with updated variables
end
%%-------------------------------------------------------------------------
%% problem dimensions
nb = length(V);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections

%% lambda selection
nlam = length(lambda) / 2;    %AAB- Number of lambdas for lamP and lamQ
lamP = lambda(1:nlam);        %AAB- From lambda vector the 1st half is for lamP
lamQ = lambda((1:nlam)+nlam); %AAB- From lambda vector the 2nd half select the lamQ

%%----- evaluate Hessian of power balance constraints -----
%% compute 2nd derivatives
%lamP
[Gp11, Gp12, Gp21, Gp22] = d2Sbus_dV2(Ybus, V, lamP, mpopt.opf.v_cartesian); 
%AAB-----------------------------------------------------------------------
[Gp15, Gp25, Gp51, Gp52, Gp55] = d2Sbus_dxBeqz2(branch, V, lamP, mpopt.opf.v_cartesian);
[Gp16, Gp26, Gp56, Gp61, Gp62, Gp65, Gp66] = d2Sbus_dxBeqv2(branch, V, lamP, mpopt.opf.v_cartesian);
[Gp17, Gp27, Gp57, Gp67, Gp71, Gp72, Gp75, Gp76, Gp77] = d2Sbus_dxsh2(branch, V, lamP, mpopt.opf.v_cartesian);
[Gp18, Gp28, Gp58, Gp68, Gp78, Gp81, Gp82, Gp85, Gp86, Gp87, Gp88] = d2Sbus_dxqtma2(branch, V, lamP, mpopt.opf.v_cartesian);
[Gp19, Gp29, Gp59, Gp69, Gp79, Gp89, Gp91, Gp92, Gp95, Gp96, Gp97, Gp98, Gp99] = d2Sbus_dxvtma2(branch, V, lamP, mpopt.opf.v_cartesian);
[Gp110, Gp210, Gp510, Gp610, Gp710, Gp810, Gp910, Gp101, Gp102, Gp105, Gp106, Gp107, Gp108, Gp109, Gp1010] = d2Sbus_dxshdp2(branch, V, lamP, mpopt.opf.v_cartesian);

%--------------------------------------------------------------------------
%lamQ
[Gq11, Gq12, Gq21, Gq22] = d2Sbus_dV2(Ybus, V, lamQ, mpopt.opf.v_cartesian);
%AAB-----------------------------------------------------------------------
[Gq15, Gq25, Gq51, Gq52, Gq55] = d2Sbus_dxBeqz2(branch, V, lamQ, mpopt.opf.v_cartesian);
[Gq16, Gq26, Gq56, Gq61, Gq62, Gq65, Gq66] = d2Sbus_dxBeqv2(branch, V, lamQ, mpopt.opf.v_cartesian);
[Gq17, Gq27, Gq57, Gq67, Gq71, Gq72, Gq75, Gq76, Gq77] = d2Sbus_dxsh2(branch, V, lamQ, mpopt.opf.v_cartesian);
[Gq18, Gq28, Gq58, Gq68, Gq78, Gq81, Gq82, Gq85, Gq86, Gq87, Gq88] = d2Sbus_dxqtma2(branch, V, lamQ, mpopt.opf.v_cartesian);
[Gq19, Gq29, Gq59, Gq69, Gq79, Gq89, Gq91, Gq92, Gq95, Gq96, Gq97, Gq98, Gq99] = d2Sbus_dxvtma2(branch, V, lamQ, mpopt.opf.v_cartesian);
[Gq110, Gq210, Gq510, Gq610, Gq710, Gq810, Gq910, Gq101, Gq102, Gq105, Gq106, Gq107, Gq108, Gq109, Gq1010] = d2Sbus_dxshdp2(branch, V, lamQ, mpopt.opf.v_cartesian);
%--------------------------------------------------------------------------

if ~mpopt.opf.v_cartesian
    %% adjust for voltage dependent loads (constant impedance part of ZIP loads)
    diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
    Sd = makeSdzip(mpc.baseMVA, mpc.bus, mpopt);
    diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
    Gp22 = Gp22 + 2 * diaglam * diagSdz;
end

%% construct Hessian
%d2G = [
%    real([Gp11 Gp12; Gp21 Gp22]) + imag([Gq11 Gq12; Gq21 Gq22]) sparse(2*nb, 2*ng);
%    sparse(2*ng, 2*nb + 2*ng)
%AAB-----------------------------------------------------------------------
Gs11=real(Gp11)+ imag(Gq11);    Gs12=real(Gp12)+ imag(Gq12);    Gs13=sparse(nb,ng);          Gs14=sparse(nb,ng);          Gs15=real(Gp15)+ imag(Gq15);          Gs16=real(Gp16)+ imag(Gq16);          Gs17=real(Gp17)+ imag(Gq17);          Gs18=real(Gp18)+ imag(Gq18);          Gs19=real(Gp19)+ imag(Gq19);          Gs110=real(Gp110)+ imag(Gq110);       
Gs21=real(Gp21)+ imag(Gq21);    Gs22=real(Gp22)+ imag(Gq22);    Gs23=sparse(nb,ng);          Gs24=sparse(nb,ng);          Gs25=real(Gp25)+ imag(Gq25);          Gs26=real(Gp26)+ imag(Gq26);          Gs27=real(Gp27)+ imag(Gq27);          Gs28=real(Gp28)+ imag(Gq28);          Gs29=real(Gp29)+ imag(Gq29);          Gs210=real(Gp210)+ imag(Gq210);
Gs31=sparse(ng,nb);             Gs32=sparse(ng,nb);             Gs33=sparse(ng,ng);          Gs34=sparse(ng,ng);          Gs35=sparse(ng,nBeqz);                Gs36=sparse(ng,nBeqv);                Gs37=sparse(ng,nPfsh);                Gs38=sparse(ng,nQtma);                Gs39=sparse(ng,nVtma);                Gs310=sparse(ng,nPfdp);
Gs41=sparse(ng,nb);             Gs42=sparse(ng,nb);             Gs43=sparse(ng,ng);          Gs44=sparse(ng,ng);          Gs45=sparse(ng,nBeqz);                Gs46=sparse(ng,nBeqv);                Gs47=sparse(ng,nPfsh);                Gs48=sparse(ng,nQtma);                Gs49=sparse(ng,nVtma);                Gs410=sparse(ng,nPfdp);
Gs51=real(Gp51)+ imag(Gq51);    Gs52=real(Gp52)+ imag(Gq52);    Gs53=sparse(nBeqz,ng);       Gs54=sparse(nBeqz,ng);       Gs55=real(Gp55)+ imag(Gq55);          Gs56=real(Gp56)+ imag(Gq56);          Gs57=real(Gp57)+ imag(Gq57);          Gs58=real(Gp58)+ imag(Gq58);          Gs59=real(Gp59)+ imag(Gq59);          Gs510=real(Gp510)+ imag(Gq510);
Gs61=real(Gp61)+ imag(Gq61);    Gs62=real(Gp62)+ imag(Gq62);    Gs63=sparse(nBeqv,ng);       Gs64=sparse(nBeqv,ng);       Gs65=real(Gp65)+ imag(Gq65);          Gs66=real(Gp66)+ imag(Gq66);          Gs67=real(Gp67)+ imag(Gq67);          Gs68=real(Gp68)+ imag(Gq68);          Gs69=real(Gp69)+ imag(Gq69);          Gs610=real(Gp610)+ imag(Gq610);
Gs71=real(Gp71)+ imag(Gq71);    Gs72=real(Gp72)+ imag(Gq72);    Gs73=sparse(nPfsh,ng);       Gs74=sparse(nPfsh,ng);       Gs75=real(Gp75)+ imag(Gq75);          Gs76=real(Gp76)+ imag(Gq76);          Gs77=real(Gp77)+ imag(Gq77);          Gs78=real(Gp78)+ imag(Gq78);          Gs79=real(Gp79)+ imag(Gq79);          Gs710=real(Gp710)+ imag(Gq710);
Gs81=real(Gp81)+ imag(Gq81);    Gs82=real(Gp82)+ imag(Gq82);    Gs83=sparse(nQtma,ng);       Gs84=sparse(nQtma,ng);       Gs85=real(Gp85)+ imag(Gq85);          Gs86=real(Gp86)+ imag(Gq86);          Gs87=real(Gp87)+ imag(Gq87);          Gs88=real(Gp88)+ imag(Gq88);          Gs89=real(Gp89)+ imag(Gq89);          Gs810=real(Gp810)+ imag(Gq810);
Gs91=real(Gp91)+ imag(Gq91);    Gs92=real(Gp92)+ imag(Gq92);    Gs93=sparse(nVtma,ng);       Gs94=sparse(nVtma,ng);       Gs95=real(Gp95)+ imag(Gq95);          Gs96=real(Gp96)+ imag(Gq96);          Gs97=real(Gp97)+ imag(Gq97);          Gs98=real(Gp98)+ imag(Gq98);          Gs99=real(Gp99)+ imag(Gq99);          Gs910=real(Gp910)+ imag(Gq910);
Gs101=real(Gp101)+ imag(Gq101); Gs102=real(Gp102)+ imag(Gq102); Gs103=sparse(nPfdp,ng);      Gs104=sparse(nPfdp,ng);      Gs105=real(Gp105)+ imag(Gq105);       Gs106=real(Gp106)+ imag(Gq106);       Gs107=real(Gp107)+ imag(Gq107);       Gs108=real(Gp108)+ imag(Gq108);       Gs109=real(Gp109)+ imag(Gq109);       Gs1010=real(Gp1010)+ imag(Gq1010);


      %Va      Vm      Pg      Qg      Beqz    Beqv    ShAng   Qtma    Vtma    Pfdp
d2G = [Gs11    Gs12    Gs13    Gs14    Gs15    Gs16    Gs17    Gs18    Gs19    Gs110;
       Gs21    Gs22    Gs23    Gs24    Gs25    Gs26    Gs27    Gs28    Gs29    Gs210;
       Gs31    Gs32    Gs33    Gs34    Gs35    Gs36    Gs37    Gs38    Gs39    Gs310;
       Gs41    Gs42    Gs43    Gs44    Gs45    Gs46    Gs47    Gs48    Gs49    Gs410;
       Gs51    Gs52    Gs53    Gs54    Gs55    Gs56    Gs57    Gs58    Gs59    Gs510;
       Gs61    Gs62    Gs63    Gs64    Gs65    Gs66    Gs67    Gs68    Gs69    Gs610;
       Gs71    Gs72    Gs73    Gs74    Gs75    Gs76    Gs77    Gs78    Gs79    Gs710;
       Gs81    Gs82    Gs83    Gs84    Gs85    Gs86    Gs87    Gs88    Gs89    Gs810;
       Gs91    Gs92    Gs93    Gs94    Gs95    Gs96    Gs97    Gs98    Gs99    Gs910;
       Gs101   Gs102   Gs103   Gs104   Gs105   Gs106   Gs107   Gs108   Gs109   Gs1010;       
       ]; %AAB Power Balance Hessian including FUBM
%--------------------------------------------------------------------------
