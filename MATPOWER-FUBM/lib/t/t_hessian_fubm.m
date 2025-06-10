function t_hessian_fubm(quiet)
%T_HESSIAN_FUBM  Numerical tests of 2nd derivative code.

%   This code compares the results from the obtained derivatives against
%   the aproximated derivatives using the finite differences method.

%   FINITE DIFFERENCES METHOD
%   This method calculates the derivatives with an aproximation as:
%   f' (x) ~~ ( f(x+h) - f(x) ) / h 
%   f''(x) ~~ ( f'(x+h) - f'(x) ) / h 

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk
%
%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(384, quiet); %AAB-initializes the global test counters (Number of Total Tests)
%casefile = 'case30';
%casefile = 'fubm_caseHVDC_qt';
%casefile = 'fubm_caseHVDC_vt';
%casefile = 'fubm_case_57_14_2MTDC_ctrls';
%casefile = 'fubm_case_30_2MTDC_ctrls_vt1_pf';
%casefile = 'fubm_case_30_2MTDC_ctrls_vt2_pf';
casefile = 'fubm_case_30_2MTDC_ctrls_vt2_pf_dp';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM- Original: idx_brch
vcart=0;
%% run powerflow to get solved case
mpopt = mpoption('verbose', 0, 'out.all', 0);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
Vr = real(V);
Vi = imag(V);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
pert = 1e-8;                                    %% perturbation factor (h) for the Finite Differences Method

%% identifier of AC/DC grids
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 | branch(:,CONV)==4) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC for zero constraint
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSCII
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
nVscL = length(iVscL); %AAB- Number of VSC with power losses

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift except VSCIII [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- find VSCIII Voltage Droop Control elements
nPfdp = length(iPfdp); %number of VSCIII Voltage Droop Control elements

%%Save the Voltage-Droop control settings though the branch (   Pf - Pfset = Kdp*(Vmf - Vmfset)  )
Kdp    = branch(:,KDP   ); %Voltage Droop Slope   setting for the branch element in p.u.
Vmfset = branch(:,VF_SET); %Voltage Droop Voltage Setting for the branch element in p.u.

%% -----  run tests for polar coordinate  -----
    %%-----  create perturbed voltages  -----
    %% polar coordinate voltages (V1=Va, V2=Vm)
        coord = 'polar';
        vv = {'aa', 'av', 'va', 'vv'};
        V1p = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
        V2p = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));

    %% -----  check d2Sbus_dV2 code  -----
    t = ' - d2Sbus_dV2 (complex power injections)';
    lam = 10 * rand(nb, 1);
    %%sparse matrices partial derivatives
    [H11, H12, H21, H22] = d2Sbus_dV2(Ybus, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    num_H11 = zeros(nb, nb);
    num_H12 = zeros(nb, nb);
    num_H21 = zeros(nb, nb);
    num_H22 = zeros(nb, nb);
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart); %dSbus_dVa
    %VaVa
    for i = 1:nb
        V1p = V;
        V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
        [dSbus_dV1_1p, dSbus_dV2_1p] = dSbus_dV(Ybus, V1p, vcart); %dSbus_dVaPertVa
        num_H11(:, i) = (dSbus_dV1_1p - dSbus_dV1).' * lam / pert; %VaVa (dSbus_dVaPertVa - dSbus_dVa)
    end
    %VaVm
    for i = 1:nb
        V2p = V;
        V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSbus_dV1_2p, dSbus_dV2_2p] = dSbus_dV(Ybus, V2p, vcart); %dSbus_dVaPertVm
        num_H12(:, i) = (dSbus_dV1_2p - dSbus_dV1).' * lam / pert; %VaVm (dSbus_dVaPertVm - dSbus_dVa)
    end
    %VmVa
    for i = 1:nb
        V1p = V;
        V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
        [dSbus_dV1_1p, dSbus_dV2_1p] = dSbus_dV(Ybus, V1p, vcart); %dSbus_dVmPertVa
        num_H21(:, i) = (dSbus_dV2_1p - dSbus_dV2).' * lam / pert; %VmVa (dSbus_dVmPertVa - dSbus_dVm)
    end
    %VmVm
    for i = 1:nb
        V2p = V;
        V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSbus_dV1_2p, dSbus_dV2_2p] = dSbus_dV(Ybus, V2p, vcart); %dSbus_dVmPertVm
        num_H22(:, i) = (dSbus_dV2_2p - dSbus_dV2).' * lam / pert; %VmVm (dSbus_dVmPertVm - dSbus_dVm)
    end
    
    t_is(full(H11), num_H11, 4, sprintf('%s - H%s%s', coord, vv{1}, t));
    t_is(full(H12), num_H12, 4, sprintf('%s - H%s%s', coord, vv{2}, t));
    
    t_is(full(H21), num_H21, 4, sprintf('%s - H%s%s', coord, vv{3}, t));
    t_is(full(H22), num_H22, 4, sprintf('%s - H%s%s', coord, vv{4}, t));


    %% -----  check d2Sbus_dxBeqz2 code  -----
    t = ' - d2Sbus_dxBeqz2 (Beqz complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G15, G25, G51, G52, G55] = d2Sbus_dxBeqz2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G15, num_G25, num_G51, num_G52, num_G55] = d2Sbus_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G15), num_G15, 4, sprintf('%s - HVaBeqz%s', coord, t));
    t_is(full(G25), num_G25, 4, sprintf('%s - HVmBeqz%s', coord, t));
    
    t_is(full(G51), num_G51, 4, sprintf('%s - HBeqzVa%s', coord, t));
    t_is(full(G52), num_G52, 4, sprintf('%s - HBeqzVm%s', coord, t));
    
    t_is(full(G55), num_G55, 4, sprintf('%s - HBeqz2 %s', coord, t));

    %% -----  check d2Sbus_dxBeqv2 code  -----
    t = ' - d2Sbus_dxBeqv2 (Beqv complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G16, G26, G56, G61, G62, G65, G66] = d2Sbus_dxBeqv2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G16, num_G26, num_G56, num_G61, num_G62, num_G65, num_G66] = d2Sbus_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G16), num_G16, 4, sprintf('%s - HVaBeqv%s'  , coord, t));
    t_is(full(G26), num_G26, 4, sprintf('%s - HVmBeqv%s'  , coord, t));
    t_is(full(G56), num_G56, 4, sprintf('%s - HBeqzBeqv%s', coord, t));
    
    t_is(full(G61), num_G61, 4, sprintf('%s - HBeqvVa%s'  , coord, t));
    t_is(full(G62), num_G62, 4, sprintf('%s - HBeqvVm%s'  , coord, t));
    t_is(full(G65), num_G65, 4, sprintf('%s - HBeqvBeqz%s', coord, t));
    
    t_is(full(G66), num_G66, 4, sprintf('%s - HBeqv2 %s'  , coord, t));

    %% -----  check d2Sbus_dxsh2 code  -----
    t = ' - d2Sbus_dxsh2 (Pf_sh complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G17, G27, G57, G67, G71, G72, G75, G76, G77] = d2Sbus_dxsh2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G17, num_G27, num_G57, num_G67, num_G71, num_G72, num_G75, num_G76, num_G77] = d2Sbus_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G17), num_G17, 4, sprintf('%s - HVaPfsh%s'  , coord, t));
    t_is(full(G27), num_G27, 4, sprintf('%s - HVmPfsh%s'  , coord, t));
    t_is(full(G57), num_G57, 4, sprintf('%s - HBeqzPfsh%s', coord, t));
    t_is(full(G67), num_G67, 4, sprintf('%s - HBeqvPfsh%s', coord, t));
    
    t_is(full(G71), num_G71, 4, sprintf('%s - HPfshVa%s'  , coord, t));
    t_is(full(G72), num_G72, 4, sprintf('%s - HPfshVm%s'  , coord, t));
    t_is(full(G75), num_G75, 4, sprintf('%s - HPfshBeqz%s', coord, t));
    t_is(full(G76), num_G76, 4, sprintf('%s - HPfshBeqv%s', coord, t));
    
    t_is(full(G77), num_G77, 4, sprintf('%s - HPfsh2 %s'  , coord, t));

    %% -----  check d2Sbus_dxqtma2 code  -----
    t = ' - d2Sbus_dxqtma2 (Qt_ma complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G18, G28, G58, G68, G78, G81, G82, G85, G86, G87, G88] = d2Sbus_dxqtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G18, num_G28, num_G58, num_G68, num_G78, num_G81, num_G82, num_G85, num_G86, num_G87, num_G88] = d2Sbus_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G18), num_G18, 4, sprintf('%s - HVaQtma%s'  , coord, t));
    t_is(full(G28), num_G28, 4, sprintf('%s - HVmQtma%s'  , coord, t));
    t_is(full(G58), num_G58, 4, sprintf('%s - HBeqzQtma%s', coord, t));
    t_is(full(G68), num_G68, 4, sprintf('%s - HBeqvQtma%s', coord, t));
    t_is(full(G78), num_G78, 4, sprintf('%s - HPfshQtma%s', coord, t));
    
    t_is(full(G81), num_G81, 4, sprintf('%s - HQtmaVa%s'  , coord, t));
    t_is(full(G82), num_G82, 4, sprintf('%s - HQtmaVm%s'  , coord, t));
    t_is(full(G85), num_G85, 4, sprintf('%s - HQtmaBeqz%s', coord, t));
    t_is(full(G86), num_G86, 4, sprintf('%s - HQtmaBeqv%s', coord, t));
    t_is(full(G87), num_G87, 4, sprintf('%s - HQtmaPfsh%s', coord, t));
    
    t_is(full(G88), num_G88, 4, sprintf('%s - HQtma2 %s'  , coord, t));
    
    %% -----  check d2Sbus_dxvtma2 code  -----
    t = ' - d2Sbus_dxvtma2 (Vt_ma complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G19, G29, G59, G69, G79, G89, G91, G92, G95, G96, G97, G98, G99] = d2Sbus_dxvtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G19, num_G29, num_G59, num_G69, num_G79, num_G89, num_G91, num_G92, num_G95, num_G96, num_G97, num_G98, num_G99] = d2Sbus_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G19), num_G19, 4, sprintf('%s - HVaVtma%s'  , coord, t));
    t_is(full(G29), num_G29, 4, sprintf('%s - HVmVtma%s'  , coord, t));
    t_is(full(G59), num_G59, 4, sprintf('%s - HBeqzVtma%s', coord, t));
    t_is(full(G69), num_G69, 4, sprintf('%s - HBeqvVtma%s', coord, t));
    t_is(full(G79), num_G79, 4, sprintf('%s - HPfshVtma%s', coord, t));
    t_is(full(G89), num_G89, 4, sprintf('%s - HQtmaVtma%s', coord, t));
    
    t_is(full(G91), num_G91, 4, sprintf('%s - HVtmaVa%s'  , coord, t));
    t_is(full(G92), num_G92, 4, sprintf('%s - HVtmaVm%s'  , coord, t));
    t_is(full(G95), num_G95, 4, sprintf('%s - HVtmaBeqz%s', coord, t));
    t_is(full(G96), num_G96, 4, sprintf('%s - HVtmaBeqv%s', coord, t));
    t_is(full(G97), num_G97, 4, sprintf('%s - HVtmaPfsh%s', coord, t));
    t_is(full(G98), num_G98, 4, sprintf('%s - HVtmaQtma%s', coord, t));
    
    t_is(full(G99), num_G99, 4, sprintf('%s - HVtma2 %s'  , coord, t));

    
    %% -----  check d2Sbus_dshdp2 code  -----
    t = ' - d2Sbus_dxshdp2 (Pf-Vf Droop Control)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G110, G210, G510, G610, G710, G810, G910, G101, G102, G105, G106, G107, G108, G109, G1010] = d2Sbus_dxshdp2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G110, num_G210, num_G510, num_G610, num_G710, num_G810, num_G910, num_G101, num_G102, num_G105, num_G106, num_G107, num_G108, num_G109, num_G1010] = d2Sbus_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G110), num_G110, 4, sprintf('%s - HVaPfdp%s'  , coord, t));
    t_is(full(G210), num_G210, 4, sprintf('%s - HVmPfdp%s'  , coord, t));
    t_is(full(G510), num_G510, 4, sprintf('%s - HBeqzPfdp%s', coord, t));
    t_is(full(G610), num_G610, 4, sprintf('%s - HBeqvPfdp%s', coord, t));
    t_is(full(G710), num_G710, 4, sprintf('%s - HPfshPfdp%s', coord, t));
    t_is(full(G810), num_G810, 4, sprintf('%s - HQtmaPfdp%s', coord, t));
    t_is(full(G910), num_G910, 4, sprintf('%s - HVtmaPfdp%s', coord, t));
    
    t_is(full(G101), num_G101, 4, sprintf('%s - HPfdpVa%s'  , coord, t));
    t_is(full(G102), num_G102, 4, sprintf('%s - HPfdpVm%s'  , coord, t));
    t_is(full(G105), num_G105, 4, sprintf('%s - HPfdpBeqz%s', coord, t));
    t_is(full(G106), num_G106, 4, sprintf('%s - HPfdpBeqv%s', coord, t));
    t_is(full(G107), num_G107, 4, sprintf('%s - HPfdpPfsh%s', coord, t));
    t_is(full(G108), num_G108, 4, sprintf('%s - HPfdpQtma%s', coord, t));
    t_is(full(G109), num_G109, 4, sprintf('%s - HPfdpVtma%s', coord, t));
    
    t_is(full(G1010), num_G1010, 4, sprintf('%s - HPfdp2 %s'  , coord, t));
    
    %% -----  check d2Sbr_dV2 code  -----
    t = ' - d2Sbr_dV2 (complex power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [Gf11, Gf12, Gf21, Gf22] = d2Sbr_dV2(Cf, Yf, V, lam, vcart);
    [Gt11, Gt12, Gt21, Gt22] = d2Sbr_dV2(Ct, Yt, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        num_Gf11(:, i) = (dSf_dV1_1p - dSf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dSf_dV2_1p - dSf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dSt_dV1_1p - dSt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dSt_dV2_1p - dSt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        num_Gf12(:, i) = (dSf_dV1_2p - dSf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dSf_dV2_2p - dSf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dSt_dV1_2p - dSt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dSt_dV2_2p - dSt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 4, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 4, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 4, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 4, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 4, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 4, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 4, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 4, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %% -----  check d2Sbr_dxBeqz2 code  -----
    t = ' - d2Sbr_dxBeqz2 (Beqz complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Sf_dxBeqz2(branch, V, lam, vcart);
    [Ht13, Ht23, Ht31, Ht32, Ht33] = d2St_dxBeqz2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33,...
     num_Ht13, num_Ht23, num_Ht31, num_Ht32, num_Ht33] = d2Sbr_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf13), num_Hf13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Hf23), num_Hf23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    
    t_is(full(Hf31), num_Hf31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Hf32), num_Hf32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    
    t_is(full(Hf33), num_Hf33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht13), num_Ht13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Ht23), num_Ht23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    
    t_is(full(Ht31), num_Ht31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Ht32), num_Ht32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    
    t_is(full(Ht33), num_Ht33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));

    %% -----  check d2Sbr_dxBeqv2 code  -----
    t = ' - d2Sbr_dxBeqv2 (Beqv complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Sf_dxBeqv2(branch, V, lam, vcart);
    [Ht14, Ht24, Ht34, Ht41, Ht42, Ht43, Ht44] = d2St_dxBeqv2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44,...
     num_Ht14, num_Ht24, num_Ht34, num_Ht41, num_Ht42, num_Ht43, num_Ht44] = d2Sbr_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf14), num_Hf14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Hf24), num_Hf24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Hf34), num_Hf34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    
    t_is(full(Hf41), num_Hf41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Hf42), num_Hf42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Hf43), num_Hf43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    
    t_is(full(Hf44), num_Hf44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht14), num_Ht14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Ht24), num_Ht24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Ht34), num_Ht34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    
    t_is(full(Ht41), num_Ht41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Ht42), num_Ht42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Ht43), num_Ht43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    
    t_is(full(Ht44), num_Ht44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    %% -----  check d2Sbr_dxPfsh2 code  -----
    t = ' - d2Sbr_dxPfsh2 (Pfsh complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = d2Sf_dxsh2(branch, V, lam, vcart);
    [Ht15, Ht25, Ht35, Ht45, Ht51, Ht52, Ht53, Ht54, Ht55] = d2St_dxsh2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55,...
     num_Ht15, num_Ht25, num_Ht35, num_Ht45, num_Ht51, num_Ht52, num_Ht53, num_Ht54, num_Ht55] = d2Sbr_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf15), num_Hf15, 4, sprintf('%s - HVaPfsh%s%s'  , coord, t, br));
    t_is(full(Hf25), num_Hf25, 4, sprintf('%s - HVmPfsh%s%s'  , coord, t, br));
    t_is(full(Hf35), num_Hf35, 4, sprintf('%s - HBeqzPfsh%s%s', coord, t, br));
    t_is(full(Hf45), num_Hf45, 4, sprintf('%s - HBeqvPfsh%s%s', coord, t, br));
    
    t_is(full(Hf51), num_Hf51, 4, sprintf('%s - HPfshVa%s%s'  , coord, t, br));
    t_is(full(Hf52), num_Hf52, 4, sprintf('%s - HPfshVm%s%s'  , coord, t, br));
    t_is(full(Hf53), num_Hf53, 4, sprintf('%s - HPfshBeqz%s%s', coord, t, br));
    t_is(full(Hf54), num_Hf54, 4, sprintf('%s - HPfshBeqv%s%s', coord, t, br));
    
    t_is(full(Hf55), num_Hf55, 4, sprintf('%s - HPfsh2 %s%s'  , coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht15), num_Ht15, 4, sprintf('%s - HVaPfsh%s%s'  , coord, t, br));
    t_is(full(Ht25), num_Ht25, 4, sprintf('%s - HVmPfsh%s%s'  , coord, t, br));
    t_is(full(Ht35), num_Ht35, 4, sprintf('%s - HBeqzPfsh%s%s', coord, t, br));
    t_is(full(Ht45), num_Ht45, 4, sprintf('%s - HBeqvPfsh%s%s', coord, t, br));
    
    t_is(full(Ht51), num_Ht51, 4, sprintf('%s - HPfshVa%s%s'  , coord, t, br));
    t_is(full(Ht52), num_Ht52, 4, sprintf('%s - HPfshVm%s%s'  , coord, t, br));
    t_is(full(Ht53), num_Ht53, 4, sprintf('%s - HPfshBeqz%s%s', coord, t, br));
    t_is(full(Ht54), num_Ht54, 4, sprintf('%s - HPfshBeqv%s%s', coord, t, br));
    
    t_is(full(Ht55), num_Ht55, 4, sprintf('%s - HPfsh2 %s%s'  , coord, t, br));

    %% -----  check d2Sbr_dxQtma2 code  -----
    t = ' - d2Sbr_dxQtma2 (Qtma complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = d2Sf_dxqtma2(branch, V, lam, vcart);
    [Ht16, Ht26, Ht36, Ht46, Ht56, Ht61, Ht62, Ht63, Ht64, Ht65, Ht66] = d2St_dxqtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66,...
     num_Ht16, num_Ht26, num_Ht36, num_Ht46, num_Ht56, num_Ht61, num_Ht62, num_Ht63, num_Ht64, num_Ht65, num_Ht66] = d2Sbr_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf16), num_Hf16, 4, sprintf('%s - HVaQtma%s%s'  , coord, t, br));
    t_is(full(Hf26), num_Hf26, 4, sprintf('%s - HVmQtma%s%s'  , coord, t, br));
    t_is(full(Hf36), num_Hf36, 4, sprintf('%s - HBeqzQtma%s%s', coord, t, br));
    t_is(full(Hf46), num_Hf46, 4, sprintf('%s - HBeqvQtma%s%s', coord, t, br));
    t_is(full(Hf56), num_Hf56, 4, sprintf('%s - HPfshQtma%s%s', coord, t, br));
    
    t_is(full(Hf61), num_Hf61, 4, sprintf('%s - HQtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf62), num_Hf62, 4, sprintf('%s - HQtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf63), num_Hf63, 4, sprintf('%s - HQtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf64), num_Hf64, 4, sprintf('%s - HQtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf65), num_Hf65, 4, sprintf('%s - HQtmaPfsh%s%s', coord, t, br));
    
    t_is(full(Hf66), num_Hf66, 4, sprintf('%s - HQtma2 %s%s'  , coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht16), num_Ht16, 4, sprintf('%s - HVaQtma%s%s'  , coord, t, br));
    t_is(full(Ht26), num_Ht26, 4, sprintf('%s - HVmQtma%s%s'  , coord, t, br));
    t_is(full(Ht36), num_Ht36, 4, sprintf('%s - HBeqzQtma%s%s', coord, t, br));
    t_is(full(Ht46), num_Ht46, 4, sprintf('%s - HBeqvQtma%s%s', coord, t, br));
    t_is(full(Ht56), num_Ht56, 4, sprintf('%s - HPfshQtma%s%s', coord, t, br));
    
    t_is(full(Ht61), num_Ht61, 4, sprintf('%s - HQtmaVa%s%s'  , coord, t, br));
    t_is(full(Ht62), num_Ht62, 4, sprintf('%s - HQtmaVm%s%s'  , coord, t, br));
    t_is(full(Ht63), num_Ht63, 4, sprintf('%s - HQtmaBeqz%s%s', coord, t, br));
    t_is(full(Ht64), num_Ht64, 4, sprintf('%s - HQtmaBeqv%s%s', coord, t, br));
    t_is(full(Ht65), num_Ht65, 4, sprintf('%s - HQtmaPfsh%s%s', coord, t, br));
    
    t_is(full(Ht66), num_Ht66, 4, sprintf('%s - HQtma2 %s%s'  , coord, t, br));
    
    %% -----  check d2Sbr_dxVtma2 code  -----
    t = ' - d2Sbr_dxVtma2 (Vtma complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = d2Sf_dxvtma2(branch, V, lam, vcart);
    [Ht17, Ht27, Ht37, Ht47, Ht57, Ht67, Ht71, Ht72, Ht73, Ht74, Ht75, Ht76, Ht77] = d2St_dxvtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf17, num_Hf27, num_Hf37, num_Hf47, num_Hf57, num_Hf67, num_Hf71, num_Hf72, num_Hf73, num_Hf74, num_Hf75, num_Hf76, num_Hf77,...
     num_Ht17, num_Ht27, num_Ht37, num_Ht47, num_Ht57, num_Ht67, num_Ht71, num_Ht72, num_Ht73, num_Ht74, num_Ht75, num_Ht76, num_Ht77] = d2Sbr_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf17), num_Hf17, 4, sprintf('%s - HVaVtma%s%s'  , coord, t, br));
    t_is(full(Hf27), num_Hf27, 4, sprintf('%s - HVmVtma%s%s'  , coord, t, br));
    t_is(full(Hf37), num_Hf37, 4, sprintf('%s - HBeqzVtma%s%s', coord, t, br));
    t_is(full(Hf47), num_Hf47, 4, sprintf('%s - HBeqvVtma%s%s', coord, t, br));
    t_is(full(Hf57), num_Hf57, 4, sprintf('%s - HPfshVtma%s%s', coord, t, br));
    t_is(full(Hf67), num_Hf67, 4, sprintf('%s - HQtmaVtma%s%s', coord, t, br));
    
    t_is(full(Hf71), num_Hf71, 4, sprintf('%s - HVtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf72), num_Hf72, 4, sprintf('%s - HVtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf73), num_Hf73, 4, sprintf('%s - HVtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf74), num_Hf74, 4, sprintf('%s - HVtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf75), num_Hf75, 4, sprintf('%s - HVtmaPfsh%s%s', coord, t, br));
    t_is(full(Hf76), num_Hf76, 4, sprintf('%s - HVtmaQtma%s%s', coord, t, br));
    
    t_is(full(Hf77), num_Hf77, 4, sprintf('%s - HVtma2 %s%s'  , coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht17), num_Ht17, 4, sprintf('%s - HVaVtma%s%s'  , coord, t, br));
    t_is(full(Ht27), num_Ht27, 4, sprintf('%s - HVmVtma%s%s'  , coord, t, br));
    t_is(full(Ht37), num_Ht37, 4, sprintf('%s - HBeqzVtma%s%s', coord, t, br));
    t_is(full(Ht47), num_Ht47, 4, sprintf('%s - HBeqvVtma%s%s', coord, t, br));
    t_is(full(Ht57), num_Ht57, 4, sprintf('%s - HPfshVtma%s%s', coord, t, br));
    t_is(full(Ht67), num_Ht67, 4, sprintf('%s - HQtmaVtma%s%s', coord, t, br));
    
    t_is(full(Ht71), num_Ht71, 4, sprintf('%s - HVtmaVa%s%s'  , coord, t, br));
    t_is(full(Ht72), num_Ht72, 4, sprintf('%s - HVtmaVm%s%s'  , coord, t, br));
    t_is(full(Ht73), num_Ht73, 4, sprintf('%s - HVtmaBeqz%s%s', coord, t, br));
    t_is(full(Ht74), num_Ht74, 4, sprintf('%s - HVtmaBeqv%s%s', coord, t, br));
    t_is(full(Ht75), num_Ht75, 4, sprintf('%s - HVtmaPfsh%s%s', coord, t, br));
    t_is(full(Ht76), num_Ht76, 4, sprintf('%s - HVtmaQtma%s%s', coord, t, br));
   
    t_is(full(Ht77), num_Ht77, 4, sprintf('%s - HVtma2 %s%s'  , coord, t, br));
    
    %% -----  check d2Sbr_dxshdp2 code  -----
    t = ' - d2Sbr_dxshdp2 (Pfdp complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = d2Sf_dxshdp2(branch, V, lam, vcart);
    [Ht18, Ht28, Ht38, Ht48, Ht58, Ht68, Ht78, Ht81, Ht82, Ht83, Ht84, Ht85, Ht86, Ht87, Ht88] = d2St_dxshdp2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf18, num_Hf28, num_Hf38, num_Hf48, num_Hf58, num_Hf68, num_Hf78, num_Hf81, num_Hf82, num_Hf83, num_Hf84, num_Hf85, num_Hf86, num_Hf87, num_Hf88,...
     num_Ht18, num_Ht28, num_Ht38, num_Ht48, num_Ht58, num_Ht68, num_Ht78, num_Ht81, num_Ht82, num_Ht83, num_Ht84, num_Ht85, num_Ht86, num_Ht87, num_Ht88] = d2Sbr_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf18), num_Hf18, 4, sprintf('%s - HVaPfdp%s%s'  , coord, t, br));
    t_is(full(Hf28), num_Hf28, 4, sprintf('%s - HVmPfdp%s%s'  , coord, t, br));
    t_is(full(Hf38), num_Hf38, 4, sprintf('%s - HBeqzPfdp%s%s', coord, t, br));
    t_is(full(Hf48), num_Hf48, 4, sprintf('%s - HBeqvPfdp%s%s', coord, t, br));
    t_is(full(Hf58), num_Hf58, 4, sprintf('%s - HPfshPfdp%s%s', coord, t, br));
    t_is(full(Hf68), num_Hf68, 4, sprintf('%s - HQtmaPfdp%s%s', coord, t, br));
    t_is(full(Hf78), num_Hf78, 4, sprintf('%s - HVtmaPfdp%s%s', coord, t, br));
    
    t_is(full(Hf81), num_Hf81, 4, sprintf('%s - HPfdpVa%s%s'  , coord, t, br));
    t_is(full(Hf82), num_Hf82, 4, sprintf('%s - HPfdpVm%s%s'  , coord, t, br));
    t_is(full(Hf83), num_Hf83, 4, sprintf('%s - HPfdpBeqz%s%s', coord, t, br));
    t_is(full(Hf84), num_Hf84, 4, sprintf('%s - HPfdpBeqv%s%s', coord, t, br));
    t_is(full(Hf85), num_Hf85, 4, sprintf('%s - HPfdpPfsh%s%s', coord, t, br));
    t_is(full(Hf86), num_Hf86, 4, sprintf('%s - HPfdpQtma%s%s', coord, t, br));
    t_is(full(Hf87), num_Hf87, 4, sprintf('%s - HPfdpVtma%s%s', coord, t, br));
    
    t_is(full(Hf88), num_Hf88, 4, sprintf('%s - HPfdp2 %s%s'  , coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht18), num_Ht18, 4, sprintf('%s - HVaPfdp%s%s'  , coord, t, br));
    t_is(full(Ht28), num_Ht28, 4, sprintf('%s - HVmPfdp%s%s'  , coord, t, br));
    t_is(full(Ht38), num_Ht38, 4, sprintf('%s - HBeqzPfdp%s%s', coord, t, br));
    t_is(full(Ht48), num_Ht48, 4, sprintf('%s - HBeqvPfdp%s%s', coord, t, br));
    t_is(full(Ht58), num_Ht58, 4, sprintf('%s - HPfshPfdp%s%s', coord, t, br));
    t_is(full(Ht68), num_Ht68, 4, sprintf('%s - HQtmaPfdp%s%s', coord, t, br));
    t_is(full(Ht78), num_Ht78, 4, sprintf('%s - HVtmaPfdp%s%s', coord, t, br));
    
    t_is(full(Ht81), num_Ht81, 4, sprintf('%s - HPfdpVa%s%s'  , coord, t, br));
    t_is(full(Ht82), num_Ht82, 4, sprintf('%s - HPfdpVm%s%s'  , coord, t, br));
    t_is(full(Ht83), num_Ht83, 4, sprintf('%s - HPfdpBeqz%s%s', coord, t, br));
    t_is(full(Ht84), num_Ht84, 4, sprintf('%s - HPfdpBeqv%s%s', coord, t, br));
    t_is(full(Ht85), num_Ht85, 4, sprintf('%s - HPfdpPfsh%s%s', coord, t, br));
    t_is(full(Ht86), num_Ht86, 4, sprintf('%s - HPfdpQtma%s%s', coord, t, br));
    t_is(full(Ht87), num_Ht87, 4, sprintf('%s - HPfdpVtma%s%s', coord, t, br));
    
    t_is(full(Ht88), num_Ht88, 4, sprintf('%s - HPfdp2 %s%s'  , coord, t, br));
    
    %% -----  check d2Abr_dV2 code  -----
    t = ' - d2Abr_dV2 (squared apparent power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [Gf11, Gf12, Gf21, Gf22] = d2Abr_dV2(d2Sf_dV2, dSf_dV1, dSf_dV2, Sf, V, lam);
    [Gt11, Gt12, Gt21, Gt22] = d2Abr_dV2(d2St_dV2, dSt_dV1, dSt_dV2, St, V, lam);
    for i = 1:nb
        V1p = V;
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va

        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        [dAf_dV1_1p, dAf_dV2_1p, dAt_dV1_1p, dAt_dV2_1p] = ...
            dAbr_dV(dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p);
        num_Gf11(:, i) = (dAf_dV1_1p - dAf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dAf_dV2_1p - dAf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dAt_dV1_1p - dAt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dAt_dV2_1p - dAt_dV2).' * lam / pert;
        
        V2p = V;
        V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        [dAf_dV1_2p, dAf_dV2_2p, dAt_dV1_2p, dAt_dV2_2p] = ...
            dAbr_dV(dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p);
        num_Gf12(:, i) = (dAf_dV1_2p - dAf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dAf_dV2_2p - dAf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dAt_dV1_2p - dAt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dAt_dV2_2p - dAt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 2.5, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 2.5, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 2.5, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 2.5, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 2.5, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 2.5, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 2.5, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 2.5, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %% -----  check d2Abr_dxfubm2 code  -----
    lam = 10 * rand(nl, 1);
    %%Annonymus functions 
    d2Sf_dBeqz2 = @(V, mu)d2Sf_dxBeqz2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqz of Sf 
    d2St_dBeqz2 = @(V, mu)d2St_dxBeqz2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqz of St 
    d2Sf_dBeqv2 = @(V, mu)d2Sf_dxBeqv2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqv of Sf 
    d2St_dBeqv2 = @(V, mu)d2St_dxBeqv2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Beqv of St 
    d2Sf_dPfsh2 = @(V, mu)d2Sf_dxsh2(branch, V, mu, vcart);     %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of Sf
    d2St_dPfsh2 = @(V, mu)d2St_dxsh2(branch, V, mu, vcart);     %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_sh of St
    d2Sf_dQtma2 = @(V, mu)d2Sf_dxqtma2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. qtma     of Sf
    d2St_dQtma2 = @(V, mu)d2St_dxqtma2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. qtma     of St
    d2Sf_dVtma2 = @(V, mu)d2Sf_dxvtma2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. vtma     of Sf
    d2St_dVtma2 = @(V, mu)d2St_dxvtma2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. vtma     of St
    d2Sf_dPfdp2 = @(V, mu)d2Sf_dxshdp2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of Sf
    d2St_dPfdp2 = @(V, mu)d2St_dxshdp2(branch, V, mu, vcart);   %AAB-Anonymus function for the 2nd derivatives w.r.t. Theta_dp of St
    
    %%First Derivatives Sbr
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1, vcart);
    [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2, vcart);
    [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1, vcart);
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, vcart);
    [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4, vcart);
    [dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3, vcart);    
    
    %%sparse matrices partial derivatives
    [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Abr_dxBeqz2(d2Sf_dBeqz2, dSf_dV1, dSf_dV2, dSf_dBeqz, Sf, V, lam);
    [Ht13, Ht23, Ht31, Ht32, Ht33] = d2Abr_dxBeqz2(d2St_dBeqz2, dSt_dV1, dSt_dV2, dSt_dBeqz, St, V, lam);
    [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Abr_dxBeqv2(d2Sf_dBeqv2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, Sf, V, lam);
    [Ht14, Ht24, Ht34, Ht41, Ht42, Ht43, Ht44] = d2Abr_dxBeqv2(d2St_dBeqv2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, St, V, lam);
    [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = d2Abr_dxsh2(d2Sf_dPfsh2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, Sf, V, lam); 
    [Ht15, Ht25, Ht35, Ht45, Ht51, Ht52, Ht53, Ht54, Ht55] = d2Abr_dxsh2(d2St_dPfsh2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, St, V, lam); 
    [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = d2Abr_dxqtma2(d2Sf_dQtma2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, Sf, V, lam); 
    [Ht16, Ht26, Ht36, Ht46, Ht56, Ht61, Ht62, Ht63, Ht64, Ht65, Ht66] = d2Abr_dxqtma2(d2St_dQtma2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, St, V, lam); 
    [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = d2Abr_dxvtma2(d2Sf_dVtma2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, dSf_dVtma, Sf, V, lam); 
    [Ht17, Ht27, Ht37, Ht47, Ht57, Ht67, Ht71, Ht72, Ht73, Ht74, Ht75, Ht76, Ht77] = d2Abr_dxvtma2(d2St_dVtma2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, dSt_dVtma, St, V, lam);   
    [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = d2Abr_dxshdp2(d2Sf_dPfdp2, dSf_dV1, dSf_dV2, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, dSf_dVtma, dSf_dPfdp, Sf, V, lam); 
    [Ht18, Ht28, Ht38, Ht48, Ht58, Ht68, Ht78, Ht81, Ht82, Ht83, Ht84, Ht85, Ht86, Ht87, Ht88] = d2Abr_dxshdp2(d2St_dPfdp2, dSt_dV1, dSt_dV2, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, dSt_dVtma, dSt_dPfdp, St, V, lam);   
    
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33,...
     num_Ht13, num_Ht23, num_Ht31, num_Ht32, num_Ht33] = d2Abr_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44,...
     num_Ht14, num_Ht24, num_Ht34, num_Ht41, num_Ht42, num_Ht43, num_Ht44] = d2Abr_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55,...
     num_Ht15, num_Ht25, num_Ht35, num_Ht45, num_Ht51, num_Ht52, num_Ht53, num_Ht54, num_Ht55] = d2Abr_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart); 
    [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66,...
     num_Ht16, num_Ht26, num_Ht36, num_Ht46, num_Ht56, num_Ht61, num_Ht62, num_Ht63, num_Ht64, num_Ht65, num_Ht66] = d2Abr_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart); 
    [num_Hf17, num_Hf27, num_Hf37, num_Hf47, num_Hf57, num_Hf67, num_Hf71, num_Hf72, num_Hf73, num_Hf74, num_Hf75, num_Hf76, num_Hf77,...
     num_Ht17, num_Ht27, num_Ht37, num_Ht47, num_Ht57, num_Ht67, num_Ht71, num_Ht72, num_Ht73, num_Ht74, num_Ht75, num_Ht76, num_Ht77] = d2Abr_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);   
    [num_Hf18, num_Hf28, num_Hf38, num_Hf48, num_Hf58, num_Hf68, num_Hf78, num_Hf81, num_Hf82, num_Hf83, num_Hf84, num_Hf85, num_Hf86, num_Hf87, num_Hf88,...
     num_Ht18, num_Ht28, num_Ht38, num_Ht48, num_Ht58, num_Ht68, num_Ht78, num_Ht81, num_Ht82, num_Ht83, num_Ht84, num_Ht85, num_Ht86, num_Ht87, num_Ht88] = d2Abr_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart);   
    
    t = ' - d2Abr_dxBeqz2 (Beqz complex power flows)';
    br = ' - "from" side';
    t_is(full(Hf13), num_Hf13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Hf23), num_Hf23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    t_is(full(Hf31), num_Hf31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Hf32), num_Hf32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    t_is(full(Hf33), num_Hf33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht13), num_Ht13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Ht23), num_Ht23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    t_is(full(Ht31), num_Ht31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Ht32), num_Ht32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    t_is(full(Ht33), num_Ht33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));
    
    t = ' - d2Abr_dxBeqv2 (Beqv complex power flows)';   
    br = ' - "from" side';
    t_is(full(Hf14), num_Hf14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Hf24), num_Hf24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Hf34), num_Hf34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    t_is(full(Hf41), num_Hf41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Hf42), num_Hf42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Hf43), num_Hf43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    t_is(full(Hf44), num_Hf44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht14), num_Ht14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Ht24), num_Ht24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Ht34), num_Ht34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    t_is(full(Ht41), num_Ht41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Ht42), num_Ht42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Ht43), num_Ht43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    t_is(full(Ht44), num_Ht44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    t = ' - d2Abr_dxPfsh2 (Pfsh complex power flows)';  
        br = ' - "from" side';
    t_is(full(Hf15), num_Hf15, 4, sprintf('%s - HVaPfsh%s%s'  , coord, t, br));
    t_is(full(Hf25), num_Hf25, 4, sprintf('%s - HVmPfsh%s%s'  , coord, t, br));
    t_is(full(Hf35), num_Hf35, 4, sprintf('%s - HBeqzPfsh%s%s', coord, t, br));
    t_is(full(Hf45), num_Hf45, 4, sprintf('%s - HBeqvPfsh%s%s', coord, t, br));
    t_is(full(Hf51), num_Hf51, 4, sprintf('%s - HPfshVa%s%s'  , coord, t, br));
    t_is(full(Hf52), num_Hf52, 4, sprintf('%s - HPfshVm%s%s'  , coord, t, br));
    t_is(full(Hf53), num_Hf53, 4, sprintf('%s - HPfshBeqz%s%s', coord, t, br));
    t_is(full(Hf54), num_Hf54, 4, sprintf('%s - HPfshBeqv%s%s', coord, t, br));
    t_is(full(Hf55), num_Hf55, 4, sprintf('%s - HPfsh2 %s%s'  , coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht15), num_Ht15, 4, sprintf('%s - HVaPfsh%s%s'  , coord, t, br));
    t_is(full(Ht25), num_Ht25, 4, sprintf('%s - HVmPfsh%s%s'  , coord, t, br));
    t_is(full(Ht35), num_Ht35, 4, sprintf('%s - HBeqzPfsh%s%s', coord, t, br));
    t_is(full(Ht45), num_Ht45, 4, sprintf('%s - HBeqvPfsh%s%s', coord, t, br));
    t_is(full(Ht51), num_Ht51, 4, sprintf('%s - HPfshVa%s%s'  , coord, t, br));
    t_is(full(Ht52), num_Ht52, 4, sprintf('%s - HPfshVm%s%s'  , coord, t, br));
    t_is(full(Ht53), num_Ht53, 4, sprintf('%s - HPfshBeqz%s%s', coord, t, br));
    t_is(full(Ht54), num_Ht54, 4, sprintf('%s - HPfshBeqv%s%s', coord, t, br));
    t_is(full(Ht55), num_Ht55, 4, sprintf('%s - HPfsh2 %s%s'  , coord, t, br));
    
    t = ' - d2Abr_dxQtma2 (Qtma complex power flows)';  
    br = ' - "from" side';
    t_is(full(Hf16), num_Hf16, 4, sprintf('%s - HVaQtma%s%s'  , coord, t, br));
    t_is(full(Hf26), num_Hf26, 4, sprintf('%s - HVmQtma%s%s'  , coord, t, br));
    t_is(full(Hf36), num_Hf36, 4, sprintf('%s - HBeqzQtma%s%s', coord, t, br));
    t_is(full(Hf46), num_Hf46, 4, sprintf('%s - HBeqvQtma%s%s', coord, t, br));
    t_is(full(Hf56), num_Hf56, 4, sprintf('%s - HPfshQtma%s%s', coord, t, br));
    t_is(full(Hf61), num_Hf61, 4, sprintf('%s - HQtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf62), num_Hf62, 4, sprintf('%s - HQtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf63), num_Hf63, 4, sprintf('%s - HQtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf64), num_Hf64, 4, sprintf('%s - HQtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf65), num_Hf65, 4, sprintf('%s - HQtmaPfsh%s%s', coord, t, br));
    t_is(full(Hf66), num_Hf66, 4, sprintf('%s - HQtma2 %s%s'  , coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht16), num_Ht16, 4, sprintf('%s - HVaQtma%s%s'  , coord, t, br));
    t_is(full(Ht26), num_Ht26, 4, sprintf('%s - HVmQtma%s%s'  , coord, t, br));
    t_is(full(Ht36), num_Ht36, 4, sprintf('%s - HBeqzQtma%s%s', coord, t, br));
    t_is(full(Ht46), num_Ht46, 4, sprintf('%s - HBeqvQtma%s%s', coord, t, br));
    t_is(full(Ht56), num_Ht56, 4, sprintf('%s - HPfshQtma%s%s', coord, t, br));
    t_is(full(Ht61), num_Ht61, 4, sprintf('%s - HQtmaVa%s%s'  , coord, t, br));
    t_is(full(Ht62), num_Ht62, 4, sprintf('%s - HQtmaVm%s%s'  , coord, t, br));
    t_is(full(Ht63), num_Ht63, 4, sprintf('%s - HQtmaBeqz%s%s', coord, t, br));
    t_is(full(Ht64), num_Ht64, 4, sprintf('%s - HQtmaBeqv%s%s', coord, t, br));
    t_is(full(Ht65), num_Ht65, 4, sprintf('%s - HQtmaPfsh%s%s', coord, t, br));
    t_is(full(Ht66), num_Ht66, 4, sprintf('%s - HQtma2 %s%s'  , coord, t, br));    
    
    t = ' - d2Abr_dxVtma2 (Vtma complex power flows)';  
    br = ' - "from" side';
    t_is(full(Hf17), num_Hf17, 4, sprintf('%s - HVaVtma%s%s'  , coord, t, br));
    t_is(full(Hf27), num_Hf27, 4, sprintf('%s - HVmVtma%s%s'  , coord, t, br));
    t_is(full(Hf37), num_Hf37, 4, sprintf('%s - HBeqzVtma%s%s', coord, t, br));
    t_is(full(Hf47), num_Hf47, 4, sprintf('%s - HBeqvVtma%s%s', coord, t, br));
    t_is(full(Hf57), num_Hf57, 4, sprintf('%s - HPfshVtma%s%s', coord, t, br));
    t_is(full(Hf67), num_Hf67, 4, sprintf('%s - HQtmaVtma%s%s', coord, t, br));
    t_is(full(Hf71), num_Hf71, 4, sprintf('%s - HVtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf72), num_Hf72, 4, sprintf('%s - HVtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf73), num_Hf73, 4, sprintf('%s - HVtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf74), num_Hf74, 4, sprintf('%s - HVtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf75), num_Hf75, 4, sprintf('%s - HVtmaPfsh%s%s', coord, t, br));
    t_is(full(Hf76), num_Hf76, 4, sprintf('%s - HVtmaQtma%s%s', coord, t, br));
    t_is(full(Hf77), num_Hf77, 4, sprintf('%s - HVtma2 %s%s'  , coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht17), num_Ht17, 4, sprintf('%s - HVaVtma%s%s'  , coord, t, br));
    t_is(full(Ht27), num_Ht27, 4, sprintf('%s - HVmVtma%s%s'  , coord, t, br));
    t_is(full(Ht37), num_Ht37, 4, sprintf('%s - HBeqzVtma%s%s', coord, t, br));
    t_is(full(Ht47), num_Ht47, 4, sprintf('%s - HBeqvVtma%s%s', coord, t, br));
    t_is(full(Ht57), num_Ht57, 4, sprintf('%s - HPfshVtma%s%s', coord, t, br));
    t_is(full(Ht67), num_Ht67, 4, sprintf('%s - HQtmaVtma%s%s', coord, t, br));
    t_is(full(Ht71), num_Ht71, 4, sprintf('%s - HVtmaVa%s%s'  , coord, t, br));
    t_is(full(Ht72), num_Ht72, 4, sprintf('%s - HVtmaVm%s%s'  , coord, t, br));
    t_is(full(Ht73), num_Ht73, 4, sprintf('%s - HVtmaBeqz%s%s', coord, t, br));
    t_is(full(Ht74), num_Ht74, 4, sprintf('%s - HVtmaBeqv%s%s', coord, t, br));
    t_is(full(Ht75), num_Ht75, 4, sprintf('%s - HVtmaPfsh%s%s', coord, t, br));
    t_is(full(Ht76), num_Ht76, 4, sprintf('%s - HVtmaQtma%s%s', coord, t, br));
    t_is(full(Ht77), num_Ht77, 4, sprintf('%s - HVtma2 %s%s'  , coord, t, br));

    t = ' - d2Abr_dxPfdp2 (Pfdp complex power flows)';  
    br = ' - "from" side';
    t_is(full(Hf18), num_Hf18, 4, sprintf('%s - HVaPfdp%s%s'  , coord, t, br));
    t_is(full(Hf28), num_Hf28, 4, sprintf('%s - HVmPfdp%s%s'  , coord, t, br));
    t_is(full(Hf38), num_Hf38, 4, sprintf('%s - HBeqzPfdp%s%s', coord, t, br));
    t_is(full(Hf48), num_Hf48, 4, sprintf('%s - HBeqvPfdp%s%s', coord, t, br));
    t_is(full(Hf58), num_Hf58, 4, sprintf('%s - HPfshPfdp%s%s', coord, t, br));
    t_is(full(Hf68), num_Hf68, 4, sprintf('%s - HQtmaPfdp%s%s', coord, t, br));
    t_is(full(Hf78), num_Hf78, 4, sprintf('%s - HVtmaPfdp%s%s', coord, t, br));    
    t_is(full(Hf81), num_Hf81, 4, sprintf('%s - HPfdpVa%s%s'  , coord, t, br));
    t_is(full(Hf82), num_Hf82, 4, sprintf('%s - HPfdpVm%s%s'  , coord, t, br));
    t_is(full(Hf83), num_Hf83, 4, sprintf('%s - HPfdpBeqz%s%s', coord, t, br));
    t_is(full(Hf84), num_Hf84, 4, sprintf('%s - HPfdpBeqv%s%s', coord, t, br));
    t_is(full(Hf85), num_Hf85, 4, sprintf('%s - HPfdpPfsh%s%s', coord, t, br));
    t_is(full(Hf86), num_Hf86, 4, sprintf('%s - HPfdpQtma%s%s', coord, t, br));
    t_is(full(Hf87), num_Hf87, 4, sprintf('%s - HPfdpVtma%s%s', coord, t, br));
    t_is(full(Hf88), num_Hf88, 4, sprintf('%s - HPfdp2 %s%s'  , coord, t, br));
    br = ' - " to " side';
    t_is(full(Ht18), num_Ht18, 4, sprintf('%s - HVaPfdp%s%s'  , coord, t, br));
    t_is(full(Ht28), num_Ht28, 4, sprintf('%s - HVmPfdp%s%s'  , coord, t, br));
    t_is(full(Ht38), num_Ht38, 4, sprintf('%s - HBeqzPfdp%s%s', coord, t, br));
    t_is(full(Ht48), num_Ht48, 4, sprintf('%s - HBeqvPfdp%s%s', coord, t, br));
    t_is(full(Ht58), num_Ht58, 4, sprintf('%s - HPfshPfdp%s%s', coord, t, br));
    t_is(full(Ht68), num_Ht68, 4, sprintf('%s - HQtmaPfdp%s%s', coord, t, br));
    t_is(full(Ht78), num_Ht78, 4, sprintf('%s - HVtmaPfdp%s%s', coord, t, br));
    t_is(full(Ht81), num_Ht81, 4, sprintf('%s - HPfdpVa%s%s'  , coord, t, br));
    t_is(full(Ht82), num_Ht82, 4, sprintf('%s - HPfdpVm%s%s'  , coord, t, br));
    t_is(full(Ht83), num_Ht83, 4, sprintf('%s - HPfdpBeqz%s%s', coord, t, br));
    t_is(full(Ht84), num_Ht84, 4, sprintf('%s - HPfdpBeqv%s%s', coord, t, br));
    t_is(full(Ht85), num_Ht85, 4, sprintf('%s - HPfdpPfsh%s%s', coord, t, br));
    t_is(full(Ht86), num_Ht86, 4, sprintf('%s - HPfdpQtma%s%s', coord, t, br));
    t_is(full(Ht87), num_Ht87, 4, sprintf('%s - HPfdpVtma%s%s', coord, t, br));
    t_is(full(Ht88), num_Ht88, 4, sprintf('%s - HPfdp2 %s%s'  , coord, t, br));
    
    %% -----  check d2Pfdp_dV2 code  -----
    t = ' - d2Pfdp_dV2 (Voltages Droop Control )';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Hf11 = zeros(nb, nb);
    num_Hf12 = zeros(nb, nb);
    num_Hf21 = zeros(nb, nb);
    num_Hf22 = zeros(nb, nb);
    
    %%Derivatives of Voltage Magnitude w.r.t. Voltage magnitude
    dVmf_dVm = sparse(zeros(nl,nb));                          % Initialize for speed [nl,nb]
    fdp = branch(iPfdp, F_BUS);                               % List of "from" buses with Voltage Droop Control [nPfdp, 1]
    Cfdp = sparse(1:nPfdp, fdp, ones(nPfdp, 1), nPfdp, nb);   % connection matrix for line & from buses with Voltage Droop Control [nPfdp, nb]
    dVmf_dVm(iPfdp,:)=Cfdp;                                   % Fill derivatives [nl, nb]
    
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    %%Partials of Pfdp w.r.t. Va 
    dPfdp_dV1 = -real(dSf_dV1);
    %%Partials of Pfdp w.r.t. Vm 
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    
    [Hf11, Hf12, Hf21, Hf22] = d2Sbr_dV2(Cf, Yf, V, lam, vcart); %Second Derivatives of Branch
    [Hf11, Hf12, Hf21, Hf22] = deal(-real(Hf11), -real(Hf12), -real(Hf21), -real(Hf22)); %Second Derivatives of Droop
    for i = 1:nb
        V1p = V;
        V2p = V;
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        dPfdp_dV1_1p = -real(dSf_dV1_1p);
        dPfdp_dV2_1p = -real(dSf_dV2_1p) + Kdp.*( dVmf_dVm ); 
        num_Hf11(:, i) = (dPfdp_dV1_1p - dPfdp_dV1).' * lam / pert;
        num_Hf21(:, i) = (dPfdp_dV2_1p - dPfdp_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        dPfdp_dV1_2p = -real(dSf_dV1_2p);
        dPfdp_dV2_2p = -real(dSf_dV2_2p) + Kdp.*( dVmf_dVm ); 
        num_Hf12(:, i) = (dPfdp_dV1_2p - dPfdp_dV1).' * lam / pert;
        num_Hf22(:, i) = (dPfdp_dV2_2p - dPfdp_dV2).' * lam / pert;
    end

    t_is(full(Hf11), num_Hf11, 4, sprintf('%s - Hf%s%s', coord, vv{1}, t));
    t_is(full(Hf12), num_Hf12, 4, sprintf('%s - Hf%s%s', coord, vv{2}, t));
    t_is(full(Hf21), num_Hf21, 4, sprintf('%s - Hf%s%s', coord, vv{3}, t));
    t_is(full(Hf22), num_Hf22, 4, sprintf('%s - Hf%s%s', coord, vv{4}, t));

    %% -----  check d2Pfdp_dxBeqz2 code  -----
    t = ' - d2Pfdp_dxBeqz2 (Beqz Droop Control)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Sf_dxBeqz2(branch, V, lam, vcart);
    [Hf13, Hf23, Hf31, Hf32, Hf33] = deal(-real(Hf13), -real(Hf23), -real(Hf31), -real(Hf32), -real(Hf33));
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33] = d2Pfdp_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf13), num_Hf13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Hf23), num_Hf23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    
    t_is(full(Hf31), num_Hf31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Hf32), num_Hf32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    
    t_is(full(Hf33), num_Hf33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));

    %% -----  check d2Pfdp_dxBeqv2 code  -----
    t = ' - d2Pfdp_dxBeqv2 (Beqv Droop Control)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Sf_dxBeqv2(branch, V, lam, vcart);
    [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = deal(-real(Hf14), -real(Hf24), -real(Hf34), -real(Hf41), -real(Hf42), -real(Hf43), -real(Hf44));
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44] = d2Pfdp_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf14), num_Hf14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Hf24), num_Hf24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Hf34), num_Hf34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    
    t_is(full(Hf41), num_Hf41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Hf42), num_Hf42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Hf43), num_Hf43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    
    t_is(full(Hf44), num_Hf44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    %% -----  check d2Sbr_dxPfsh2 code  -----
    t = ' - d2Pfdp_dxPfsh2 (Pfsh Droop Control)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = d2Sf_dxsh2(branch, V, lam, vcart);
    [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = deal(-real(Hf15), -real(Hf25), -real(Hf35), -real(Hf45), -real(Hf51), -real(Hf52), -real(Hf53), -real(Hf54), -real(Hf55));
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf15, num_Hf25, num_Hf35, num_Hf45, num_Hf51, num_Hf52, num_Hf53, num_Hf54, num_Hf55] = d2Pfdp_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf15), num_Hf15, 4, sprintf('%s - HVaPfsh%s%s'  , coord, t, br));
    t_is(full(Hf25), num_Hf25, 4, sprintf('%s - HVmPfsh%s%s'  , coord, t, br));
    t_is(full(Hf35), num_Hf35, 4, sprintf('%s - HBeqzPfsh%s%s', coord, t, br));
    t_is(full(Hf45), num_Hf45, 4, sprintf('%s - HBeqvPfsh%s%s', coord, t, br));
    
    t_is(full(Hf51), num_Hf51, 4, sprintf('%s - HPfshVa%s%s'  , coord, t, br));
    t_is(full(Hf52), num_Hf52, 4, sprintf('%s - HPfshVm%s%s'  , coord, t, br));
    t_is(full(Hf53), num_Hf53, 4, sprintf('%s - HPfshBeqz%s%s', coord, t, br));
    t_is(full(Hf54), num_Hf54, 4, sprintf('%s - HPfshBeqv%s%s', coord, t, br));
    
    t_is(full(Hf55), num_Hf55, 4, sprintf('%s - HPfsh2 %s%s'  , coord, t, br));
    
    %% -----  check d2Pfdp_dxQtma2 code  -----
    t = ' - d2Pfdp_dxQtma2 (Qtma Control Droop)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = d2Sf_dxqtma2(branch, V, lam, vcart);
    [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = deal(-real(Hf16), -real(Hf26), -real(Hf36), -real(Hf46), -real(Hf56), -real(Hf61), -real(Hf62), -real(Hf63), -real(Hf64), -real(Hf65), -real(Hf66));
  
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf16, num_Hf26, num_Hf36, num_Hf46, num_Hf56, num_Hf61, num_Hf62, num_Hf63, num_Hf64, num_Hf65, num_Hf66] = d2Pfdp_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf16), num_Hf16, 4, sprintf('%s - HVaQtma%s%s'  , coord, t, br));
    t_is(full(Hf26), num_Hf26, 4, sprintf('%s - HVmQtma%s%s'  , coord, t, br));
    t_is(full(Hf36), num_Hf36, 4, sprintf('%s - HBeqzQtma%s%s', coord, t, br));
    t_is(full(Hf46), num_Hf46, 4, sprintf('%s - HBeqvQtma%s%s', coord, t, br));
    t_is(full(Hf56), num_Hf56, 4, sprintf('%s - HPfshQtma%s%s', coord, t, br));
    
    t_is(full(Hf61), num_Hf61, 4, sprintf('%s - HQtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf62), num_Hf62, 4, sprintf('%s - HQtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf63), num_Hf63, 4, sprintf('%s - HQtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf64), num_Hf64, 4, sprintf('%s - HQtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf65), num_Hf65, 4, sprintf('%s - HQtmaPfsh%s%s', coord, t, br));
    
    t_is(full(Hf66), num_Hf66, 4, sprintf('%s - HQtma2 %s%s'  , coord, t, br));
    
    %% -----  check d2Pfdp_dxVtma2 code  -----
    t = ' - d2Pfdp_dxVtma2 (Vtma Droop Control)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = d2Sf_dxvtma2(branch, V, lam, vcart);
    [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = deal(-real(Hf17), -real(Hf27), -real(Hf37), -real(Hf47), -real(Hf57), -real(Hf67), -real(Hf71), -real(Hf72), -real(Hf73), -real(Hf74), -real(Hf75), -real(Hf76), -real(Hf77));
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf17, num_Hf27, num_Hf37, num_Hf47, num_Hf57, num_Hf67, num_Hf71, num_Hf72, num_Hf73, num_Hf74, num_Hf75, num_Hf76, num_Hf77] = d2Pfdp_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf17), num_Hf17, 4, sprintf('%s - HVaVtma%s%s'  , coord, t, br));
    t_is(full(Hf27), num_Hf27, 4, sprintf('%s - HVmVtma%s%s'  , coord, t, br));
    t_is(full(Hf37), num_Hf37, 4, sprintf('%s - HBeqzVtma%s%s', coord, t, br));
    t_is(full(Hf47), num_Hf47, 4, sprintf('%s - HBeqvVtma%s%s', coord, t, br));
    t_is(full(Hf57), num_Hf57, 4, sprintf('%s - HPfshVtma%s%s', coord, t, br));
    t_is(full(Hf67), num_Hf67, 4, sprintf('%s - HQtmaVtma%s%s', coord, t, br));
    
    t_is(full(Hf71), num_Hf71, 4, sprintf('%s - HVtmaVa%s%s'  , coord, t, br));
    t_is(full(Hf72), num_Hf72, 4, sprintf('%s - HVtmaVm%s%s'  , coord, t, br));
    t_is(full(Hf73), num_Hf73, 4, sprintf('%s - HVtmaBeqz%s%s', coord, t, br));
    t_is(full(Hf74), num_Hf74, 4, sprintf('%s - HVtmaBeqv%s%s', coord, t, br));
    t_is(full(Hf75), num_Hf75, 4, sprintf('%s - HVtmaPfsh%s%s', coord, t, br));
    t_is(full(Hf76), num_Hf76, 4, sprintf('%s - HVtmaQtma%s%s', coord, t, br));
    
    t_is(full(Hf77), num_Hf77, 4, sprintf('%s - HVtma2 %s%s'  , coord, t, br));
    
    %% -----  check d2Pfdp_dxshdp2 code  -----
    t = ' - d2Pfdp_dxshdp2 (Pfdp Droop Control)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = d2Sf_dxshdp2(branch, V, lam, vcart);
    [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = deal(-real(Hf18), -real(Hf28), -real(Hf38), -real(Hf48), -real(Hf58), -real(Hf68), -real(Hf78), -real(Hf81), -real(Hf82), -real(Hf83), -real(Hf84), -real(Hf85), -real(Hf86), -real(Hf87), -real(Hf88));

    %%compute numerically to compare (Finite Differences Method)
    [num_Hf18, num_Hf28, num_Hf38, num_Hf48, num_Hf58, num_Hf68, num_Hf78, num_Hf81, num_Hf82, num_Hf83, num_Hf84, num_Hf85, num_Hf86, num_Hf87, num_Hf88] = d2Pfdp_dxshdp2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf18), num_Hf18, 4, sprintf('%s - HVaPfdp%s%s'  , coord, t, br));
    t_is(full(Hf28), num_Hf28, 4, sprintf('%s - HVmPfdp%s%s'  , coord, t, br));
    t_is(full(Hf38), num_Hf38, 4, sprintf('%s - HBeqzPfdp%s%s', coord, t, br));
    t_is(full(Hf48), num_Hf48, 4, sprintf('%s - HBeqvPfdp%s%s', coord, t, br));
    t_is(full(Hf58), num_Hf58, 4, sprintf('%s - HPfshPfdp%s%s', coord, t, br));
    t_is(full(Hf68), num_Hf68, 4, sprintf('%s - HQtmaPfdp%s%s', coord, t, br));
    t_is(full(Hf78), num_Hf78, 4, sprintf('%s - HVtmaPfdp%s%s', coord, t, br));
    
    t_is(full(Hf81), num_Hf81, 4, sprintf('%s - HPfdpVa%s%s'  , coord, t, br));
    t_is(full(Hf82), num_Hf82, 4, sprintf('%s - HPfdpVm%s%s'  , coord, t, br));
    t_is(full(Hf83), num_Hf83, 4, sprintf('%s - HPfdpBeqz%s%s', coord, t, br));
    t_is(full(Hf84), num_Hf84, 4, sprintf('%s - HPfdpBeqv%s%s', coord, t, br));
    t_is(full(Hf85), num_Hf85, 4, sprintf('%s - HPfdpPfsh%s%s', coord, t, br));
    t_is(full(Hf86), num_Hf86, 4, sprintf('%s - HPfdpQtma%s%s', coord, t, br));
    t_is(full(Hf87), num_Hf87, 4, sprintf('%s - HPfdpVtma%s%s', coord, t, br));
    
    t_is(full(Hf88), num_Hf88, 4, sprintf('%s - HPfdp2 %s%s'  , coord, t, br));
    
t_end;
