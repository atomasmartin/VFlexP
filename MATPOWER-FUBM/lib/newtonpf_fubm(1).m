function [V, branch, converged, i] = newtonpf_fubm(branch, Sbus, getYbus, V0, ref, pv, pq, baseMVA, mpopt)
%NEWTONPF_FUBM  Solves the AC/DC power flow using a full Newton's method.
%   [V, BRANCH, CONVERGED, I] = NEWTONPF_FUBM(BRANCH, SBUS, GETYBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves AC/DC power flow using a full Newton-Raphson method, given the
%   following inputs:
%       BRANCH  - branch matrix (extended version for FUBM)
%       SBUS    - handle to function that returns the complex bus power
%                 injection vector (for all buses), given the bus voltage
%                 magnitude vector (for all buses)
%       GETYBUS - handle to function that returns the admitance matrix Ybus
%       V0      - initial vector of complex bus voltages
%       REF     - bus index of reference bus (voltage ang reference & gen slack)
%       PV      - vector of bus indices for PV buses
%       PQ      - vector of bus indices for PQ buses
%       BASEMVA - Power base in MVA
%       MPOPT   - (optional) MATPOWER option struct, used to set the
%                 termination tolerance, maximum number of iterations, and
%                 output options (see MPOPTION for details).
%
%   The bus voltage vector contains the set point for generator
%   (including ref bus) buses, and the reference angle of the swing
%   bus, as well as an initial guess for remaining magnitudes and
%   angles.
%
%   Returns the final complex voltages, the updated branch matrix, a flag 
%   which indicates whether it converged or not, and the number of 
%   iterations performed.
%
%   See also RUNPF, NEWTONPF.

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%% define named indices branch matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM-extra fields for FUBM- Original: idx_brch

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol        = mpopt.pf.tol;
max_it     = mpopt.pf.nr.max_it; %FUBM- Maximum number of iterations
if max_it < 25
    max_it = 25; %NR with FUBM controls or AC/DC grids sometimes takes more than 4 iterations. 
end
lin_solver  = mpopt.pf.nr.lin_solver;
%FUBM-----------------------------------------------------------------------
if mpopt.verbose > 1
    conv_graph = 1; %FUBM- Plot convergence graph option
else
    conv_graph = 0; %FUBM- Do not Plot convergence graph option
end
%% Identify FUBM formulation
if (size(branch,2) < KDP) 
    fubm = 0; %Its not a fubm formulation
    error('newtonpf_fubm: There is missing data in the branch matrix to run Power Flows using FUBM formulation. Add the data or try "newtonpf"')
else
    fubm = 1; %Its a fubm formulation
end 
%--------------------------------------------------------------------------

%% Identify elements
%%FUBM----------------------------------------------------------------------
%%Identify power control elements (Converters and transformers included)
iBeqz = find( (branch(:,CONV)==1 | branch(:,CONV)==3 ) & (branch(:,BR_STATUS)~=0) ); %FUBM- Location of the branch elements with Qf control by Beq to meet the zero constraint. (Converters type I for zero constraint control) Qf = 0
iPfsh = find( (branch(:,PF  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, SH_MIN )~=-360 | branch(:, SH_MAX )~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4) ); %FUBM- Location of the branch elements with Pf control by theta_shift to meet the setting. (Converters and Phase Shifter Transformers, but no VSCIII)
iPtsh = find( (branch(:,PT  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, SH_MIN )~=-360 | branch(:, SH_MAX )~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4) ); %FUBM- Location of the branch elements with Pt control by theta_shift to meet the setting. (Converters and Phase Shifter Transformers, but no VSCIII)
iQfma = find( (branch(:,QF  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VF_SET)==0 & branch(:,CONV)==0 ); %FUBM- Location of the branch elements with Qf control by ma to meet the setting. (Transformers) avoid all VSC
iQtma = find( (branch(:,QT  )~=0) & (branch(:,BR_STATUS)~=0) & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 & branch(:,QF  )==0 ); %FUBM- Location of the branch elements with Qt control by ma to meet the setting. (Converters and Transformers)           %FUBM- Location of the branch elements with Qt control by ma to meet the setting. (Converters and Transformers)
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %FUBM- Find VSC with active PWM Losses Calculation [nVscL,1]
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)

%Identify voltage control elements and their buses (Converters and transformers included)
iBeqv = find( (branch(:,CONV)==2) & (branch(:,BR_STATUS)~=0) & (branch(:,VF_SET)~=0) ); %FUBM- Location of the branch elements with Vf control by Beq to meet the setting. (Converters type II for Vdc control) Vf = Vfset
iVtma = find((branch(:,BR_STATUS)~=0) & (branch(:,VT_SET)~=0) & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) ); %FUBM- Location of the branch elements with Vt control by ma  to meet the setting. (Converters and Transformers)
VfBeqbus = branch(iBeqv, F_BUS); %FUBM- Saves the "from" bus identifier for Vf controlled by Beq (Converters type II for Vdc control)
Vtmabus  = branch(iVtma , T_BUS); %FUBM- Saves the "to"   bus identifier for Vt controlled by ma  (Converters and Transformers)

%%Save the Power to control through the branch
Pfset = branch(:,PF)/baseMVA; %Power constraint for the branch element in p.u.
Qfset = branch(:,QF)/baseMVA; %Power constraint for the branch element in p.u.
Qtset = branch(:,QT)/baseMVA; %Power constraint for the branch element in p.u.

%%Save the Voltage-Droop control settings though the branch (   Pf - Pfset = Kdp*(Vmf - Vmfset)  )
Kdp    = branch(:,KDP   ); %Voltage Droop Slope   setting for the branch element in p.u.
Vmfset = branch(:,VF_SET); %Voltage Droop Voltage Setting for the branch element in p.u.
%%-------------------------------------------------------------------------

%% initialize V
converged = 0; %FUBM- Convergence flag
i = 0;         %FUBM- iterations
V = V0;        %FUBM- V  initial condition
Va = angle(V); %FUBM- Va initial condition
Vm = abs(V);   %FUBM- Vm initial condition

%% set up indexing for updating state variables
npv = length(pv);
npq = length(pq);

%% variable dimensions
%%FUBM-----------------------------------------------------------------------
nPfsh = length(iPfsh); %FUBM- number of Pf controlled elements by theta_shift
nQfma = length(iQfma); %FUBM- number of Qf controlled elements by ma
nQtma = length(iQtma); %FUBM- number of Qt controlled elements by ma
nVtma = length(iVtma); %FUBM- number of Vt controlled elements by ma
nBeqz = length(iBeqz); %FUBM- number of Qf controlled elements by Beq
nBeqv = length(iBeqv); %FUBM- number of Vf controlled elements by Beq
nVscL = length(iVscL); %FUBM- Number of VSC with active PWM Losses Calculation  
nPfdp = length(iPfdp); %FUBM- Number of VSC with Voltage Droop Control by theta_shift

nVfBeqbus = length(VfBeqbus); %FUBM- number of buses for Vf controlled by Beq 
nVtmabus  = length(Vtmabus ); %FUBM- number of buses for Vt controlled by ma  
%--------------------------------------------------------------------------

%% variables dimensions in Jacobian
%%FUBM----------------------------------------------------------------------
j1  = 1;          j2   = npv;              %% j1 :j2  - V angle of pv buses (bus)
j3  = j2  + 1;    j4   = j2  + npq;        %% j3 :j4  - V angle of pq buses (bus)
j5  = j4  + 1;    j6   = j4  + npq;        %% j5 :j6  - V mag   of pq buses (bus)
j7  = j6  + 1;    j8   = j6  + nPfsh;      %% j7 :j8  - ShiftAngle of VSC and PST (branch)
j9  = j8  + 1;    j10  = j8  + nQfma;      %% j9 :j10 - ma of Qf Controlled Transformers (branch)
jA1 = j10 + 1;    jA2  = j10 + nBeqz;      %% j11:j12 - Beq of VSC for Zero Constraint (branch)
jA3 = jA2 + 1;    jA4  = jA2 + nVfBeqbus;  %% j13:j14 - Beq of VSC for Vdc  Constraint (bus)
jA5 = jA4 + 1;    jA6  = jA4 + nVtmabus;   %% j15:j16 - ma of VSC and Transformers for Vt Control (bus)
jA7 = jA6 + 1;    jA8  = jA6 + nQtma;      %% j17:j18 - ma of VSC and Transformers for Qt Control (branch)
jA9 = jA8 + 1;    jB0  = jA8 + nPfdp;      %% j19:j20 - ShiftAngle of VSC for Qt Control (branch)
%%-------------------------------------------------------------------------

%% Create addmitance matrix to initialise
%%FUBM----------------------------------------------------------------------
[Ybus, Yf, Yt] = getYbus(branch); %FUBM- Obtains Ybus, Yf, Yt from branch
%%-------------------------------------------------------------------------
%% compute branch power flows
%br = find(branch(:, BR_STATUS));  %%FUBM- in-service branches
brf=branch(:, F_BUS);              %%FUBM- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
brt=branch(:, T_BUS);              %%FUBM- to   bus index of all the branches, brt=branch(br, T_BUS); %For in-service branches
If= Yf(:, :) * V;                  %%FUBM- complex current injected at "from" bus, Yf(br, :) * V; For in-service branches 
It= Yt(:, :) * V;                  %%FUBM- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
Sf = V(brf) .* conj(If);           %%FUBM- complex power injected at "from" bus
St = V(brt) .* conj(It);           %%FUBM- complex power injected at "to"   bus
%% compute VSC Power Loss
if nVscL
    PLoss_IEC = branch(iVscL,ALPH3).*((abs(It(iVscL))).^2) + branch(iVscL,ALPH2).*((abs(It(iVscL)))) + branch(iVscL,ALPH1); %%FUBM- Standard IEC 62751-2 Ploss Correction for VSC losses 
    branch(iVscL,GSW) = PLoss_IEC./(abs(V(brf(iVscL))).^2);    %%FUBM- VSC Gsw
end
%% evaluate F(x0)
%%FUBM----------------------------------------------------------------------
mis = V .* conj(Ybus * V) - Sbus(Vm);     %FUBM- F1(x0) & F2(x0) Power balance mismatch
misPfsh = real(Sf(iPfsh)) - Pfset(iPfsh); %FUBM- F3(x0) Pf control mismatch  
misQfma = imag(Sf(iQfma)) - Qfset(iQfma); %FUBM- F4(x0) Qf control mismatch 
misBeqz = imag(Sf(iBeqz)) - 0;            %FUBM- F5(x0) Qf control mismatch 
misBeqv = imag(mis(VfBeqbus));            %FUBM- F6(x0) Vf control mismatch 
misVtma = imag(mis(Vtmabus ));            %FUBM- F7(x0) Vt control mismatch 
misQtma = imag(St(iQtma)) - Qtset(iQtma); %FUBM- F8(x0) Qt control mismatch 
misPfdp = -real(Sf(iPfdp)) + Pfset(iPfdp) + Kdp(iPfdp).* (  Vm(brf(iPfdp)) - Vmfset(iPfdp) );%FUBM- F9(x0) Pf control mismatch, Droop Pf - Pfset = Kdp*(Vmf - Vmfset)
%%-------------------------------------------------------------------------

%% Create F vector
%%FUBM----------------------------------------------------------------------
F = [   real(mis([pv; pq])); %FUBM- F1(x0) Power balance mismatch - Va
        imag(mis(pq));       %FUBM- F2(x0) Power balance mismatch - Vm
        misPfsh;             %FUBM- F3(x0) Pf control    mismatch - Theta_shift
        misQfma;             %FUBM- F4(x0) Qf control    mismatch - ma
        misBeqz;             %FUBM- F5(x0) Qf control    mismatch - Beq
        misBeqv;             %FUBM- F6(x0) Vf control    mismatch - Beq
        misVtma;             %FUBM- F7(x0) Vt control    mismatch - ma
        misQtma;             %FUBM- F8(x0) Qt control    mismatch - ma
        misPfdp;             %FUBM- F9(x0) Pf control    mismatch - Theta_shift Droop        
        ];
%%-------------------------------------------------------------------------    

%% check tolerance
normF = norm(F, inf); %FUBM- Calculates the mismatch error of the vector F, norm = max( abs(F) ) 
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = length(F);
    if nx <= 10 || have_fcn('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% do Newton iterations
while (~converged && i < max_it) %Repeat until convergence or the maximum iterations limit is reached
    %% update iteration counter
    i = i + 1;
    
    %% Convergence graph
    %%FUBM------------------------------------------------------------------
    if conv_graph
        %Saves the value of the variable to be ploted
        fmis(i,:)=F';                      %FUBM- Saves iterative values (mismatch)
        if nQtma 
            fQtma (i,:)=branch(iQtma,TAP); %FUBM- Saves iterative values (tap or VSC modulation amplitude)
        end
        if nVtma 
            fVtma (i,:)=branch(iVtma,TAP); %FUBM- Saves iterative values (tap or VSC modulation amplitude)
        end
    end
    %%---------------------------------------------------------------------
    %% evaluate Jacobian
    %FUBM-------------------------------------------------------------------
    %F1(x), F2(x), F6(x) and F7(x) Partial derivatives Power Balance w.r.t. x
    [dSbus_dVa, dSbus_dVm, dSbus_dPfsh, dSbus_dQfma,dSbus_dBeqz,...
        dSbus_dBeqv, dSbus_dVtma, dSbus_dQtma, dSbus_dPfdp] = dSbus_dx(Ybus, branch, V, 0);
    [dummy, neg_dSd_dVm] = Sbus(Vm);
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    %F3(x), F4(x), F5(x) and F8(x) Partial derivatives Sf and St w.r.t. x
    [dSf_dVa, dSf_dVm, dSf_dPfsh, dSf_dQfma, dSf_dBeqz,...
        dSf_dBeqv, dSf_dVtma, dSf_dQtma, dSf_dPfdp,...
        dSt_dVa, dSt_dVm, dSt_dPfsh, dSt_dQfma, dSt_dBeqz,...
        dSt_dBeqv, dSt_dVtma, dSt_dQtma, dSt_dPfdp] = dSbr_dx(branch, Yf, Yt, V, 0);
    %F9(x) Partial derivatives Droop Control w.r.t. x
    [dPfdp_dVa, dPfdp_dVm, dPfdp_dPfsh, dPfdp_dQfma, dPfdp_dBeqz,...
        dPfdp_dBeqv, dPfdp_dVtma, dPfdp_dQtma, dPfdp_dPfdp] = dPfdp_dx(branch, Yf, Yt, V, 1, 0); %The 1 indicates that is for Power Flows and not for OPF
       
    %----------------------------------------------------------------------
    
    %% fill Jacobian
    %FUBM-------------------------------------------------------------------
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));      %avoid Slack
    j12 = real(dSbus_dVm([pv; pq], pq));            %avoid Slack
    j13 = real(dSbus_dPfsh([pv; pq],:));            %avoid Slack
    j14 = real(dSbus_dQfma([pv; pq],:));            %avoid Slack
    j15 = real(dSbus_dBeqz([pv; pq],:));            %avoid Slack
    j16 = real(dSbus_dBeqv([pv; pq],:));            %avoid Slack
    j17 = real(dSbus_dVtma([pv; pq],:));            %avoid Slack
    j18 = real(dSbus_dQtma([pv; pq],:));            %avoid Slack
    j19 = real(dSbus_dPfdp([pv; pq],:));            %avoid Slack
    
    j21 = imag(dSbus_dVa([pq], [pv; pq]));          %avoid Slack and pv
    j22 = imag(dSbus_dVm([pq], pq));                %avoid Slack and pv
    j23 = imag(dSbus_dPfsh([pq],:));                %avoid Slack and pv
    j24 = imag(dSbus_dQfma([pq],:));                %avoid Slack and pv
    j25 = imag(dSbus_dBeqz([pq],:));                %avoid Slack and pv
    j26 = imag(dSbus_dBeqv([pq],:));                %avoid Slack and pv
    j27 = imag(dSbus_dVtma([pq],:));                %avoid Slack and pv
    j28 = imag(dSbus_dQtma([pq],:));                %avoid Slack and pv
    j29 = imag(dSbus_dPfdp([pq],:));                %avoid Slack and pv
    
    j31 = real(dSf_dVa(iPfsh,[pv; pq]));            %Only Pf control elements iPfsh
    j32 = real(dSf_dVm(iPfsh,pq));                  %Only Pf control elements iPfsh
    j33 = real(dSf_dPfsh(iPfsh,:));                 %Only Pf control elements iPfsh
    j34 = real(dSf_dQfma(iPfsh,:));                 %Only Pf control elements iPfsh
    j35 = real(dSf_dBeqz(iPfsh,:));                 %Only Pf control elements iPfsh
    j36 = real(dSf_dBeqv(iPfsh,:));                 %Only Pf control elements iPfsh
    j37 = real(dSf_dVtma(iPfsh,:));                 %Only Pf control elements iPfsh
    j38 = real(dSf_dQtma(iPfsh,:));                 %Only Pf control elements iPfsh
    j39 = real(dSf_dPfdp(iPfsh,:));                 %Only Pf control elements iPfsh
       
    j41 = imag(dSf_dVa(iQfma,[pv; pq]));            %Only Qf control elements iQfma
    j42 = imag(dSf_dVm(iQfma,pq));                  %Only Qf control elements iQfma
    j43 = imag(dSf_dPfsh(iQfma,:));                 %Only Qf control elements iQfma
    j44 = imag(dSf_dQfma(iQfma,:));                 %Only Qf control elements iQfma
    j45 = imag(dSf_dBeqz(iQfma,:));                 %Only Qf control elements iQfma
    j46 = imag(dSf_dBeqv(iQfma,:));                 %Only Qf control elements iQfma
    j47 = imag(dSf_dVtma(iQfma,:));                 %Only Qf control elements iQfma
    j48 = imag(dSf_dQtma(iQfma,:));                 %Only Qf control elements iQfma
    j49 = imag(dSf_dPfdp(iQfma,:));                 %Only Qf control elements iQfma
        
    j51 = imag(dSf_dVa(iBeqz,[pv; pq]));            %Only Qf control elements iQfbeq
    j52 = imag(dSf_dVm(iBeqz,pq));                  %Only Qf control elements iQfbeq
    j53 = imag(dSf_dPfsh(iBeqz,:));                 %Only Qf control elements iQfbeq
    j54 = imag(dSf_dQfma(iBeqz,:));                 %Only Qf control elements iQfbeq
    j55 = imag(dSf_dBeqz(iBeqz,:));                 %Only Qf control elements iQfbeq
    j56 = imag(dSf_dBeqv(iBeqz,:));                 %Only Qf control elements iQfbeq
    j57 = imag(dSf_dVtma(iBeqz,:));                 %Only Qf control elements iQfbeq
    j58 = imag(dSf_dQtma(iBeqz,:));                 %Only Qf control elements iQfbeq
    j59 = imag(dSf_dPfdp(iBeqz,:));                 %Only Qf control elements iQfbeq
    
    j61 = imag(dSbus_dVa([VfBeqbus],[pv; pq]));  %Only Vf control elements iVfbeq
    j62 = imag(dSbus_dVm([VfBeqbus],pq));        %Only Vf control elements iVfbeq
    j63 = imag(dSbus_dPfsh([VfBeqbus],:));       %Only Vf control elements iVfbeq
    j64 = imag(dSbus_dQfma([VfBeqbus],:));       %Only Vf control elements iVfbeq
    j65 = imag(dSbus_dBeqz([VfBeqbus],:));       %Only Vf control elements iVfbeq
    j66 = imag(dSbus_dBeqv([VfBeqbus],:));       %Only Vf control elements iVfbeq 
    j67 = imag(dSbus_dVtma([VfBeqbus],:));       %Only Vf control elements iVfbeq
    j68 = imag(dSbus_dQtma([VfBeqbus],:));       %Only Vf control elements iVfbeq
    j69 = imag(dSbus_dPfdp([VfBeqbus],:));       %Only Vf control elements iVfbeq
    
    j71 = imag(dSbus_dVa([Vtmabus],[pv; pq]));   %Only Vt control elements iVtma
    j72 = imag(dSbus_dVm([Vtmabus],pq));         %Only Vt control elements iVtma
    j73 = imag(dSbus_dPfsh([Vtmabus],:));        %Only Vt control elements iVtma
    j74 = imag(dSbus_dQfma([Vtmabus],:));        %Only Vt control elements iVtma 
    j75 = imag(dSbus_dBeqz([Vtmabus],:));        %Only Vt control elements iVtma
    j76 = imag(dSbus_dBeqv([Vtmabus],:));        %Only Vt control elements iVtma  
    j77 = imag(dSbus_dVtma([Vtmabus],:));        %Only Vt control elements iVtma
    j78 = imag(dSbus_dQtma([Vtmabus],:));        %Only Vt control elements iVtma
    j79 = imag(dSbus_dPfdp([Vtmabus],:));        %Only Vt control elements iVtma
    
    j81 = imag(dSt_dVa(iQtma,[pv; pq]));           %Only Qt control elements iQtma
    j82 = imag(dSt_dVm(iQtma,pq));                 %Only Qt control elements iQtma
    j83 = imag(dSt_dPfsh(iQtma,:));                %Only Qt control elements iQtma
    j84 = imag(dSt_dQfma(iQtma,:));                %Only Qt control elements iQtma
    j85 = imag(dSt_dBeqz(iQtma,:));                %Only Qt control elements iQtma
    j86 = imag(dSt_dBeqv(iQtma,:));                %Only Qt control elements iQtma
    j87 = imag(dSt_dVtma(iQtma,:));                %Only Qt control elements iQtma
    j88 = imag(dSt_dQtma(iQtma,:));                %Only Qt control elements iQtma
    j89 = imag(dSt_dPfdp(iQtma,:));                %Only Qt control elements iQtma

    j91 =      dPfdp_dVa(iPfdp,[pv; pq]);         %Only Droop control elements iPfdp
    j92 =      dPfdp_dVm(iPfdp,pq);               %Only Droop control elements iPfdp
    j93 =      dPfdp_dPfsh(iPfdp,:);              %Only Droop control elements iPfdp
    j94 =      dPfdp_dQfma(iPfdp,:);              %Only Droop control elements iPfdp
    j95 =      dPfdp_dBeqz(iPfdp,:);              %Only Droop control elements iPfdp
    j96 =      dPfdp_dBeqv(iPfdp,:);              %Only Droop control elements iPfdp
    j97 =      dPfdp_dVtma(iPfdp,:);              %Only Droop control elements iPfdp
    j98 =      dPfdp_dQtma(iPfdp,:);              %Only Droop control elements iPfdp
    j99 =      dPfdp_dPfdp(iPfdp,:);              %Only Droop control elements iPfdp
    
    %Jacobian
    J = [   j11 j12 j13 j14 j15 j16 j17 j18 j19;
            j21 j22 j23 j24 j25 j26 j27 j28 j29;   
            j31 j32 j33 j34 j35 j36 j37 j38 j39;
            j41 j42 j43 j44 j45 j46 j47 j48 j49;
            j51 j52 j53 j54 j55 j56 j57 j58 j59;
            j61 j62 j63 j64 j65 j66 j67 j68 j69;
            j71 j72 j73 j74 j75 j76 j77 j78 j79;
            j81 j82 j83 j84 j85 j86 j87 j88 j89;
            j91 j92 j93 j94 j95 j96 j97 j98 j99;
            ]; %FUBM-Jacobian Matrix
    %%---------------------------------------------------------------------
    %% compute update step
   dx = mplinsolve(J, -F, lin_solver);
   %dx = -( J \ F ); %Matpower 6
    %% update state variables
    %%FUBM------------------------------------------------------------------
    %voltage Va, Vm
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va); %%FUBM- Reconstruct V
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    %Theta_shift  PST, VSCI, VSCII
    if nPfsh 
        branch(iPfsh,SHIFT) = branch(iPfsh,SHIFT) + (dx(j7:j8).*(180/pi));
    end
    %%Qfma (tap)
    if nQfma 
        branch(iQfma,TAP) = branch(iQfma,TAP) + (dx(j9:j10));
        for k=1:nQfma
            if  branch(iQfma(k),TAP) > branch(iQfma(k),TAP_MAX)
                branch(iQfma(k),TAP) = branch(iQfma(k),TAP_MAX); %FUBM- Fix maximum tap
            elseif branch(iQfma(k),TAP) < branch(iQfma(k),TAP_MIN)
                branch(iQfma(k),TAP) = branch(iQfma(k),TAP_MIN); %FUBM- Fix minimum tap
            end
        end
    end
    %%Beqz
    if nBeqz 
        branch(iBeqz,BEQ)= branch(iBeqz,BEQ) + dx(jA1:jA2);
    end
    %%Beqv
    if nBeqv 
        branch(iBeqv,BEQ)= branch(iBeqv,BEQ) + dx(jA3:jA4);
    end
    %%Vtma (tap)
    if nVtma 
        branch(iVtma,TAP) = branch(iVtma,TAP) + (dx(jA5:jA6));
        for k=1:nVtma
            if  branch(iVtma(k),TAP) > branch(iVtma(k),TAP_MAX)
                branch(iVtma(k),TAP) = branch(iVtma(k),TAP_MAX); %FUBM- Fix maximum tap
            elseif branch(iVtma(k),TAP) < branch(iVtma(k),TAP_MIN)
                branch(iVtma(k),TAP) = branch(iVtma(k),TAP_MIN); %FUBM- Fix minimum tap
            end
        end
    end
    %%Qtma (tap)
    if nQtma 
        branch(iQtma,TAP)= branch(iQtma,TAP) + (dx(jA7:jA8));
        for k=1:nQtma
            if branch(iQtma(k),TAP) > branch(iQtma(k),TAP_MAX)
                branch(iQtma(k),TAP) = branch(iQtma(k),TAP_MAX);
            elseif branch(iQtma(k),TAP) < branch(iQtma(k),TAP_MIN)
                branch(iQtma(k),TAP) = branch(iQtma(k),TAP_MIN);
            end
        end
    end
    %Theta_shift VSCIII Droop Control
    if nPfdp 
        branch(iPfdp,SHIFT) = branch(iPfdp,SHIFT) + (dx(jA9:jB0).*(180/pi));
    end
    %%---------------------------------------------------------------------
    %% recalculate Ybus with updated variables
    %%FUBM------------------------------------------------------------------
    [Ybus, Yf, Yt] = getYbus(branch); %FUBM- Obtains Ybus, Yf, Yt
    %%---------------------------------------------------------------------
%% compute branch power flows
%br = find(branch(:, BR_STATUS));  %%FUBM- in-service branches
brf=branch(:, F_BUS);              %%FUBM- from bus index of all the branches, brf=branch(br, F_BUS); %For in-service branches 
brt=branch(:, T_BUS);              %%FUBM- to   bus index of all the branches, brt=branch(br, T_BUS); %For in-service branches
If= Yf(:, :) * V;                  %%FUBM- complex current injected at "from" bus, Yf(br, :) * V; For in-service branches 
It= Yt(:, :) * V;                  %%FUBM- complex current injected at "to"   bus, Yt(br, :) * V; For in-service branches 
Sf = V(brf) .* conj(If);           %%FUBM- complex power injected at "from" bus
St = V(brt) .* conj(It);           %%FUBM- complex power injected at "to"   bus
%% compute VSC Power Loss
if nVscL
    PLoss_IEC = branch(iVscL,ALPH3).*((abs(It(iVscL))).^2) + branch(iVscL,ALPH2).*((abs(It(iVscL)))) + branch(iVscL,ALPH1); %%FUBM- Standard IEC 62751-2 Ploss Correction for VSC losses 
    branch(iVscL,GSW) = PLoss_IEC./(abs(V(brf(iVscL))).^2);    %%FUBM- VSC Gsw
end
%% evaluate F(x0)
%%FUBM----------------------------------------------------------------------
mis = V .* conj(Ybus * V) - Sbus(Vm);     %FUBM- F1(x0) & F2(x0) Power balance mismatch
misPfsh = real(Sf(iPfsh)) - Pfset(iPfsh); %FUBM- F3(x0) Pf control mismatch  
misQfma = imag(Sf(iQfma)) - Qfset(iQfma); %FUBM- F4(x0) Qf control mismatch 
misBeqz = imag(Sf(iBeqz)) - 0;            %FUBM- F5(x0) Qf control mismatch 
misBeqv = imag(mis(VfBeqbus));            %FUBM- F6(x0) Vf control mismatch 
misVtma = imag(mis(Vtmabus ));            %FUBM- F7(x0) Vt control mismatch 
misQtma = imag(St(iQtma)) - Qtset(iQtma); %FUBM- F8(x0) Qt control mismatch 
misPfdp = -real(Sf(iPfdp)) + Pfset(iPfdp) + Kdp(iPfdp).* (  Vm(brf(iPfdp)) - Vmfset(iPfdp) );%FUBM- F9(x0) Pf control mismatch, Droop Pf - Pfset = Kdp*(Vmf - Vmfset)
%%-------------------------------------------------------------------------

%% Create F vector
%%FUBM----------------------------------------------------------------------
F = [   real(mis([pv; pq])); %FUBM- F1(x0) Power balance mismatch - Va
        imag(mis(pq));       %FUBM- F2(x0) Power balance mismatch - Vm
        misPfsh;             %FUBM- F3(x0) Pf control    mismatch - Theta_shift
        misQfma;             %FUBM- F4(x0) Qf control    mismatch - ma
        misBeqz;             %FUBM- F5(x0) Qf control    mismatch - Beq
        misBeqv;             %FUBM- F6(x0) Vf control    mismatch - Beq
        misVtma;             %FUBM- F7(x0) Vt control    mismatch - ma
        misQtma;             %FUBM- F8(x0) Qt control    mismatch - ma
        misPfdp;             %FUBM- F9(x0) Pf control    mismatch - Theta_shift Droop
        ];
%%------------------------------------------------------------------------- 
    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nNewton''s method power flow converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow did not converge in %d iterations.\n', i);
    end
end
%% Plot graphs
if conv_graph
    %Saves the last value of the variable to be ploted
    fmis(i+1,:)=F';                      %FUBM- Saves iterative values (mismatch)
    if nQtma 
        fQtma (i+1,:)=branch(iQtma,TAP); %FUBM- Saves iterative values (tap or VSC modulation amplitude)
    end
    if nVtma 
        fVtma (i+1,:)=branch(iVtma,TAP); %FUBM- Saves iterative values (tap or VSC modulation amplitude)
    end
    %Ploting
    Graphsfmis(i, fmis);
    if iQtma
        GraphsQtma(i, fQtma); 
    end
    if iVtma
        GraphsVtma(i, fVtma); 
    end
end

