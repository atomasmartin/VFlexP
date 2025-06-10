function [results] = sim_fubm
%% SIM_FUBM Quick Simulation of some cases for the FUBM 
clear all;
clc; 
%Options
% mpopt = mpoption('opf.ac.solver', 'KNITRO','knitro.tol_x',1e-10,'knitro.tol_f',1e-4,'opf.violation',1e-6,'opf.start',0);
% mpopt = mpoption('opf.ac.solver', 'MIPS', 'mips.max_it',5000,'opf.violation',1e-6,'opf.start',0);
 mpopt = mpoption('opf.ac.solver', 'FMINCON','opf.violation',1e-6,'opf.start',0);
% mpopt = mpoption('opf.ac.solver', 'IPOPT','opf.violation',1e-6,'opf.start',0);
mpopt = mpoption(mpopt, 'verbose', 2);

%Run OPF
%[results] = runopf('fubm_caseHVDC_qt',mpopt);
%[results] = runopf('fubm_caseHVDC_vt',mpopt);
%[results] = runopf('fubm_case_57_14_2MTDC_ctrls',mpopt);
%[results] = runopf('fubm_case_30_2MTDC_ctrls_vt1_pf',mpopt);
%[results] = runopf('fubm_case_30_2MTDC_ctrls_vt2_pf',mpopt);
%[results] = runopf('fubm_case_30_2MTDC_ctrls_vt2_pf_dp',mpopt);       
[results] = runopf('fubm_case1354pegase_2MTDC_ctrls_pf_qt_dp',mpopt); %Only Power Flow

%Constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM-extra fields for FUBM

%Locating controls and AC/DC grids
iBeqz = find ((results.branch(:,CONV)==1 | results.branch(:,CONV)==3 | results.branch(:,CONV)==4) & results.branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC for Zero Constraint control size[nBeqz,1]
iBeqv = find (results.branch(:,CONV)==2 & results.branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC for Zero Constraint control size[nBeqz,1]
iPfsh = find (results.branch(:,PF)  ~=0 & results.branch(:, BR_STATUS)==1 & (results.branch(:, SH_MIN)~=-360 | results.branch(:, SH_MAX)~=360)& (results.branch(:, CONV)~=3) & (results.branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
iQtma = find (results.branch(:,QT)  ~=0 & results.branch(:, BR_STATUS)==1 & results.branch(:, VT_SET)==0 & (results.branch(:, TAP_MIN)~= results.branch(:, TAP_MAX)) ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
iVtma = find (results.branch(:, BR_STATUS)==1 & results.branch(:, VT_SET)~=0 & (results.branch(:, TAP_MIN)~= results.branch(:, TAP_MAX)) ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
iVscL = find (results.branch(:,CONV)~=0 & results.branch(:, BR_STATUS)==1 & (results.branch(:, ALPH1)~=0 | results.branch(:, ALPH2)~=0 | results.branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
iPfdp = find( (results.branch(:,VF_SET)~=0) & (results.branch(:,KDP)~=0) & (results.branch(:, BR_STATUS)~=0) & (results.branch(:, SH_MIN)~=-360 | results.branch(:, SH_MAX)~=360) & (results.branch(:, CONV)==3 | results.branch(:, CONV)==4) );

%Display final Results
disp('iBeqz')
disp('----FROM----TO----------BEQZ------QF-----')
disp(results.branch(iBeqz,[F_BUS, T_BUS, BEQ, QF]))
disp('iBeqv')
disp('----FROM----TO----------BEQV----')
disp(results.branch(iBeqv,[F_BUS, T_BUS, BEQ]))
disp('iPfsh')
disp('----FROM----TO----------SHANG------PF-----')
disp(results.branch(iPfsh,[F_BUS, T_BUS, SHIFT, PF]))
disp('iQtma')
disp('----FROM----TO------------MA-------QT-----')
disp(results.branch(iQtma,[F_BUS, T_BUS, TAP, QT]))
disp('iVtma')
disp('----FROM----TO------------MA----')
disp(results.branch(iVtma,[F_BUS, T_BUS, TAP]))
disp('iVscL')
disp('----FROM----TO------------GSW----')
disp(results.branch(iVscL,[F_BUS, T_BUS, GSW]))
disp('iPfdp')
disp('----FROM----TO----------SHANG------PF-----')
disp(results.branch(iPfdp,[F_BUS, T_BUS, SHIFT, PF]))
