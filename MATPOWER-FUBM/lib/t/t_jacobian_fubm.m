function t_jacobian_fubm(quiet)
%T_JACOBIAN_FUBM  Numerical tests of partial derivative code.

%   This code compares the results from the obtained derivatives against
%   the aproximated derivatives using the finite differences method.

%   FINITE DIFFERENCES METHOD
%   This method calculates the derivatives with an aproximation as:
%   f'(x) ~~ ( f(x+h) - f(x) ) / h 

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

t_begin(54, quiet); %AAB-initializes the global test counters (Number of Total Tests)

%casefile = 'fubm_caseHVDC_qt';
%casefile = 'fubm_caseHVDC_vt';
%casefile = 'fubm_case_57_14_2MTDC_ctrls';
%casefile = 'fubm_case_30_2MTDC_ctrls_vt1_pf';
casefile = 'fubm_case_30_2MTDC_ctrls_vt2_pf_dp';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM- Original: idx_brch

%% run powerflow to get solved case
mpopt = mpoption('verbose', 0, 'out.all', 0);
mpc = loadcase(casefile);
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.  Sbus = Sbusg - Sbusd;
Ybus_full   = full(Ybus);
Yf_full     = full(Yf);
Yt_full     = full(Yt);
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
Vr = real(V);
Vi = imag(V);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);             %                                            [1,2,...,nb]
VV = V * ones(1, nb);       %% Voltages are repeated in each column VV = [V,V,..., V]; size [nb,nb]
SS = Sbus * ones(1, nb);    %% Sbus      is repeated in each column SS = [Sbus,Sbus,...,Sbus]; size [nb,nb] where: Sbus = Sbusg - Sbusd; 
pert = 1e-8;                %% perturbation factor (h) for the Finite Differences Method

%%Save the Voltage-Droop control settings though the branch (   Pf - Pfset = Kdp*(Vmf - Vmfset)  )
Kdp    = branch(:,KDP   ); %Voltage Droop Slope   setting for the branch element in p.u.
Vmfset = branch(:,VF_SET); %Voltage Droop Voltage Setting for the branch element in p.u.


%%-----  run tests for polar coordinates -----
    coord = 'polar';
    vv = {'a', 'm'}; %a for angle, m for magnitude
    V1p = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb)))); %Perturbed V for Va+pert
    V2p = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb)); %Perturbed V for Vm+pert

    %% -----  check dSbus_dx code  -----
    %%sparse matrices derivatives
    [dSbus_dV1, dSbus_dV2, dSbus_dPfsh, dSbus_dQfma,dSbus_dBeqz,...
        dSbus_dBeqv, dSbus_dVtma, dSbus_dQtma, dSbus_dPfdp] = dSbus_dx(Ybus, branch, V, 0);
    dSbus_dV1_sp   = full(dSbus_dV1);
    dSbus_dV2_sp   = full(dSbus_dV2);
    dSbus_dPfsh_sp = full(dSbus_dPfsh);
    dSbus_dQfma_sp = full(dSbus_dQfma);
    dSbus_dBeqz_sp = full(dSbus_dBeqz);
    dSbus_dBeqv_sp = full(dSbus_dBeqv);
    dSbus_dVtma_sp = full(dSbus_dVtma);
    dSbus_dQtma_sp = full(dSbus_dQtma);
    dSbus_dPfdp_sp = full(dSbus_dPfdp);
    
    %%compute numerically to compare (Finite Differences Method)
    %Voltages
    num_dSbus_dV1 = full( (V1p .* conj(Ybus * V1p) - VV .* conj(Ybus * VV)) / pert );
    num_dSbus_dV2 = full( (V2p .* conj(Ybus * V2p) - VV .* conj(Ybus * VV)) / pert );
    %FUBM extra variables
    [num_dSbus_dPfsh, num_dSbus_dQfma,num_dSbus_dBeqz,...
        num_dSbus_dBeqv, num_dSbus_dVtma, num_dSbus_dQtma, num_dSbus_dPfdp] = dSbus_dxPert(baseMVA, bus, branch, V, pert, 0); %Size of each derivatives [nb, nXXxx]
    
    num_dSbus_dPfsh_sp = full(num_dSbus_dPfsh);
    num_dSbus_dQfma_sp = full(num_dSbus_dQfma);
    num_dSbus_dBeqz_sp = full(num_dSbus_dBeqz);
    num_dSbus_dBeqv_sp = full(num_dSbus_dBeqv);
    num_dSbus_dVtma_sp = full(num_dSbus_dVtma);
    num_dSbus_dQtma_sp = full(num_dSbus_dQtma);
    num_dSbus_dPfdp_sp = full(num_dSbus_dPfdp);
    %Results comparison
    t_is(dSbus_dV1_sp, num_dSbus_dV1, 5, sprintf('%s - dSbus_dV%s (sparse)', coord, vv{1}));
    t_is(dSbus_dV2_sp, num_dSbus_dV2, 5, sprintf('%s - dSbus_dV%s (sparse)', coord, vv{2}));
    
    t_is(dSbus_dPfsh_sp, num_dSbus_dPfsh_sp, 5, sprintf('%s - dSbus_dPfsh (sparse)', coord)); 
    t_is(dSbus_dQfma_sp, num_dSbus_dQfma_sp, 5, sprintf('%s - dSbus_dQfma (sparse)', coord));
    t_is(dSbus_dBeqz_sp, num_dSbus_dBeqz_sp, 5, sprintf('%s - dSbus_dBeqz (sparse)', coord));
    t_is(dSbus_dBeqv_sp, num_dSbus_dBeqv_sp, 5, sprintf('%s - dSbus_dBeqv (sparse)', coord));    
    t_is(dSbus_dVtma_sp, num_dSbus_dVtma_sp, 5, sprintf('%s - dSbus_dVtma (sparse)', coord)); 
    t_is(dSbus_dQtma_sp, num_dSbus_dQtma_sp, 5, sprintf('%s - dSbus_dQtma (sparse)', coord)); 
    t_is(dSbus_dPfdp_sp, num_dSbus_dPfdp_sp, 5, sprintf('%s - dSbus_dPfdp (sparse)', coord));     
    %% -----  check dSbr_dx code  -----
    %%sparse matrices
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, 0);
    [dSf_dV1, dSf_dV2, dSf_dPfsh, dSf_dQfma, dSf_dBeqz,...
        dSf_dBeqv, dSf_dVtma, dSf_dQtma, dSf_dPfdp,...
        dSt_dV1, dSt_dV2, dSt_dPfsh, dSt_dQfma, dSt_dBeqz,...
        dSt_dBeqv, dSt_dVtma, dSt_dQtma, dSt_dPfdp] = dSbr_dx(branch, Yf, Yt, V, 0);
    
    dSf_dV1_sp   = full(dSf_dV1);    dSt_dV1_sp   = full(dSt_dV1);
    dSf_dV2_sp   = full(dSf_dV2);    dSt_dV2_sp   = full(dSt_dV2);
    dSf_dPfsh_sp = full(dSf_dPfsh);  dSt_dPfsh_sp = full(dSt_dPfsh);
    dSf_dQfma_sp = full(dSf_dQfma);  dSt_dQfma_sp = full(dSt_dQfma);
    dSf_dBeqz_sp = full(dSf_dBeqz);  dSt_dBeqz_sp = full(dSt_dBeqz);    
    dSf_dBeqv_sp = full(dSf_dBeqv);  dSt_dBeqv_sp = full(dSt_dBeqv);  
    dSf_dVtma_sp = full(dSf_dVtma);  dSt_dVtma_sp = full(dSt_dVtma);
    dSf_dQtma_sp = full(dSf_dQtma);  dSt_dQtma_sp = full(dSt_dQtma);
    dSf_dPfdp_sp = full(dSf_dPfdp);  dSt_dPfdp_sp = full(dSt_dPfdp);
    
    %%compute numerically to compare (Finite Differences Method)
    %Voltages "from" and "to"
    V1pf = V1p(f,:);
    V2pf = V2p(f,:);
    V1pt = V1p(t,:);
    V2pt = V2p(t,:);
    
    %Obtain f(x) and f(x+h) for Voltages
    Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * VV);
    St2 = (V(t)*ones(1,nb)) .* conj(Yt * VV);
    S1pf = V1pf .* conj(Yf * V1p);
    S2pf = V2pf .* conj(Yf * V2p);
    S1pt = V1pt .* conj(Yt * V1p);
    S2pt = V2pt .* conj(Yt * V2p);

    %Voltages
    num_dPfdp_dV1 = full( (S1pf - Sf2) / pert );
    num_dPfdp_dV2 = full( (S2pf - Sf2) / pert );
    num_dSt_dV1 = full( (S1pt - St2) / pert );
    num_dSt_dV2 = full( (S2pt - St2) / pert );
    
    %FUBM extra variables
    [num_dSf_dPfsh, num_dSf_dQfma, num_dSf_dBeqz,...
     num_dSf_dBeqv, num_dSf_dVtma, num_dSf_dQtma, num_dSf_dPfdp,...
     num_dSt_dPfsh, num_dSt_dQfma, num_dSt_dBeqz,...
     num_dSt_dBeqv, num_dSt_dVtma, num_dSt_dQtma, num_dSt_dPfdp] = dSbr_dxPert(baseMVA, bus, branch, Yf, Yt, V, pert, 0); %Size of each derivatives [nb, nXXxx]
    
    num_dSf_dPfsh_sp = full(num_dSf_dPfsh);     num_dSt_dPfsh_sp = full(num_dSt_dPfsh);
    num_dSf_dQfma_sp = full(num_dSf_dQfma);     num_dSt_dQfma_sp = full(num_dSt_dQfma);
    num_dSf_dBeqz_sp = full(num_dSf_dBeqz);     num_dSt_dBeqz_sp = full(num_dSt_dBeqz);
    num_dSf_dBeqv_sp = full(num_dSf_dBeqv);     num_dSt_dBeqv_sp = full(num_dSt_dBeqv);
    num_dSf_dVtma_sp = full(num_dSf_dVtma);     num_dSt_dVtma_sp = full(num_dSt_dVtma);
    num_dSf_dQtma_sp = full(num_dSf_dQtma);     num_dSt_dQtma_sp = full(num_dSt_dQtma);
    num_dSf_dPfdp_sp = full(num_dSf_dPfdp);     num_dSt_dPfdp_sp = full(num_dSt_dPfdp);
    %Results Comparison
    t_is(dSf_dV1_sp, num_dPfdp_dV1, 5, sprintf('%s - dSf_dV%s (sparse)', coord, vv{1}));
    t_is(dSf_dV2_sp, num_dPfdp_dV2, 5, sprintf('%s - dSf_dV%s (sparse)', coord, vv{2}));
    t_is(dSf_dPfsh_sp, num_dSf_dPfsh_sp, 5, sprintf('%s - dSf_dPfsh (sparse)', coord));
    t_is(dSf_dQfma_sp, num_dSf_dQfma_sp, 5, sprintf('%s - dSf_dQfma (sparse)', coord));
    t_is(dSf_dBeqz_sp, num_dSf_dBeqz_sp, 5, sprintf('%s - dSf_dBeqz (sparse)', coord));
    t_is(dSf_dBeqv_sp, num_dSf_dBeqv_sp, 5, sprintf('%s - dSf_dBeqv (sparse)', coord));
    t_is(dSf_dVtma_sp, num_dSf_dVtma_sp, 5, sprintf('%s - dSf_dVtma (sparse)', coord));
    t_is(dSf_dQtma_sp, num_dSf_dQtma_sp, 5, sprintf('%s - dSf_dQtma (sparse)', coord));       
    t_is(dSf_dPfdp_sp, num_dSf_dPfdp_sp, 5, sprintf('%s - dSf_dPfdp (sparse)', coord));
    
    t_is(dSt_dV1_sp, num_dSt_dV1, 5, sprintf('%s - dSt_dV%s (sparse)', coord, vv{1}));
    t_is(dSt_dV2_sp, num_dSt_dV2, 5, sprintf('%s - dSt_dV%s (sparse)', coord, vv{2}));
    t_is(dSt_dPfsh_sp, num_dSt_dPfsh_sp, 5, sprintf('%s - dSt_dPfsh (sparse)', coord));
    t_is(dSt_dQfma_sp, num_dSt_dQfma_sp, 5, sprintf('%s - dSt_dQfma (sparse)', coord));
    t_is(dSt_dBeqz_sp, num_dSt_dBeqz_sp, 5, sprintf('%s - dSt_dBeqz (sparse)', coord));
    t_is(dSt_dBeqv_sp, num_dSt_dBeqv_sp, 5, sprintf('%s - dSt_dBeqv (sparse)', coord));
    t_is(dSt_dVtma_sp, num_dSt_dVtma_sp, 5, sprintf('%s - dSt_dVtma (sparse)', coord));
    t_is(dSt_dQtma_sp, num_dSt_dQtma_sp, 5, sprintf('%s - dSt_dQtma (sparse)', coord)); 
    t_is(dSt_dPfdp_sp, num_dSt_dPfdp_sp, 5, sprintf('%s - dSt_dPfdp (sparse)', coord));
    
    %% -----  check dAbr_dV code  -----
    %%sparse matrices
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    dAf_dV1_sp = full(dAf_dV1);
    dAf_dV2_sp = full(dAf_dV2);
    dAt_dV1_sp = full(dAt_dV1);
    dAt_dV2_sp = full(dAt_dV2);

    %%compute numerically to compare
    num_dAf_dV1 = full( (abs(S1pf).^2 - abs(Sf2).^2) / pert );
    num_dAf_dV2 = full( (abs(S2pf).^2 - abs(Sf2).^2) / pert );
    num_dAt_dV1 = full( (abs(S1pt).^2 - abs(St2).^2) / pert );
    num_dAt_dV2 = full( (abs(S2pt).^2 - abs(St2).^2) / pert );

    t_is(dAf_dV1_sp, num_dAf_dV1, 4, sprintf('%s - dAf_dV%s (sparse)', coord, vv{1}));
    t_is(dAf_dV2_sp, num_dAf_dV2, 4, sprintf('%s - dAf_dV%s (sparse)', coord, vv{2}));
    t_is(dAt_dV1_sp, num_dAt_dV1, 4, sprintf('%s - dAt_dV%s (sparse)', coord, vv{1}));
    t_is(dAt_dV2_sp, num_dAt_dV2, 4, sprintf('%s - dAt_dV%s (sparse)', coord, vv{2}));

    %% -----  check dAbr_dsh code  -----
    %%sparse matrices
    [dAf_dPfsh, dAt_dPfsh] = ...
                        dAbr_dsh(dSf_dPfsh, dSt_dPfsh, Sf, St);
    [dAf_dPfdp, dAt_dPfdp] = ...
                        dAbr_dsh(dSf_dPfdp, dSt_dPfdp, Sf, St);
    dAf_dPfsh_sp = full(dAf_dPfsh);
    dAt_dPfsh_sp = full(dAt_dPfsh);
    dAf_dPfdp_sp = full(dAf_dPfdp);
    dAt_dPfdp_sp = full(dAt_dPfdp);
    
    %% compute numerically to compare
    [num_dAf_dPfsh, num_dAt_dPfsh] = ...
                    dAbr_dshPert(baseMVA, bus, branch, V, 1, pert, 0);
    [num_dAf_dPfdp, num_dAt_dPfdp] = ...
                    dAbr_dshPert(baseMVA, bus, branch, V, 3, pert, 0);
    
    t_is(dAf_dPfsh_sp, num_dAf_dPfsh, 4, sprintf('%s - dAf_dPfsh (sparse)', coord));
    t_is(dAt_dPfsh_sp, num_dAt_dPfsh, 4, sprintf('%s - dAt_dPfsh (sparse)', coord));
    t_is(dAf_dPfdp_sp, num_dAf_dPfdp, 4, sprintf('%s - dAf_dPfdp (sparse)', coord));
    t_is(dAt_dPfdp_sp, num_dAt_dPfdp, 4, sprintf('%s - dAt_dPfdp (sparse)', coord));
    
        %% -----  check dAbr_dBeqx code  -----
    %%sparse matrices
    [dAf_dBeqz, dAt_dBeqz] = ...
                        dAbr_dBeq(dSf_dBeqz, dSt_dBeqz, Sf, St);
    [dAf_dBeqv, dAt_dBeqv] = ...
                        dAbr_dBeq(dSf_dBeqv, dSt_dBeqv, Sf, St);                    
    dAf_dBeqz_sp = full(dAf_dBeqz);
    dAt_dBeqz_sp = full(dAt_dBeqz);
    dAf_dBeqv_sp = full(dAf_dBeqv);
    dAt_dBeqv_sp = full(dAt_dBeqv);

    %% compute numerically to compare
    [num_dAf_dBeqz, num_dAt_dBeqz] = ...
                    dAbr_dBeqPert(baseMVA, bus, branch, V, 1, pert, 0);
    [num_dAf_dBeqv, num_dAt_dBeqv] = ...
                    dAbr_dBeqPert(baseMVA, bus, branch, V, 2, pert, 0);
                
    t_is(dAf_dBeqz_sp, num_dAf_dBeqz, 4, sprintf('%s - dAf_dBeqz (sparse)', coord));
    t_is(dAt_dBeqz_sp, num_dAt_dBeqz, 4, sprintf('%s - dAt_dBeqz (sparse)', coord));
    t_is(dAf_dBeqv_sp, num_dAf_dBeqv, 4, sprintf('%s - dAf_dBeqv (sparse)', coord));
    t_is(dAt_dBeqv_sp, num_dAt_dBeqv, 4, sprintf('%s - dAt_dBeqv (sparse)', coord));  
    
    %% -----  check dAbr_dma code  -----
    %%sparse matrices
    [dAf_dQfma, dAt_dQfma] = ...
                        dAbr_dma(dSf_dQfma, dSt_dQfma, Sf, St);    
    [dAf_dQtma, dAt_dQtma] = ...
                        dAbr_dma(dSf_dQtma, dSt_dQtma, Sf, St);
    [dAf_dVtma, dAt_dVtma] = ...
                        dAbr_dma(dSf_dVtma, dSt_dVtma, Sf, St);                  
    dAf_dQfma_sp = full(dAf_dQfma);
    dAt_dQfma_sp = full(dAt_dQfma);
    dAf_dQtma_sp = full(dAf_dQtma);
    dAt_dQtma_sp = full(dAt_dQtma);
    dAf_dVtma_sp = full(dAf_dVtma);
    dAt_dVtma_sp = full(dAt_dVtma);

    %% compute numerically to compare
    [num_dAf_dQfma, num_dAt_dQfma] = ...
                    dAbr_dmaPert(baseMVA, bus, branch, V, 1, pert, 0);    
    [num_dAf_dQtma, num_dAt_dQtma] = ...
                    dAbr_dmaPert(baseMVA, bus, branch, V, 2, pert, 0);          
    [num_dAf_dVtma, num_dAt_dVtma] = ...
                    dAbr_dmaPert(baseMVA, bus, branch, V, 4, pert, 0);  

    t_is(dAf_dQfma_sp, num_dAf_dQfma, 4, sprintf('%s - dAf_dQfma (sparse)', coord));
    t_is(dAt_dQfma_sp, num_dAt_dQfma, 4, sprintf('%s - dAt_dQfma (sparse)', coord));
    t_is(dAf_dQtma_sp, num_dAf_dQtma, 4, sprintf('%s - dAf_dQtma (sparse)', coord));
    t_is(dAt_dQtma_sp, num_dAt_dQtma, 4, sprintf('%s - dAt_dQtma (sparse)', coord));
    t_is(dAf_dVtma_sp, num_dAf_dVtma, 4, sprintf('%s - dAf_dVtma (sparse)', coord));
    t_is(dAt_dVtma_sp, num_dAt_dVtma, 4, sprintf('%s - dAt_dVtma (sparse)', coord));
    
    %% -----  check dPfdp_dx code  -----
    %       Pfdp     = -Pf     + Pfset + Kdp.*(Vmf     - Vmfset)
    %       PfdpPert = -PfPert + Pfset + Kdp.*(VmfPert - Vmfset)    
    %       dPfdp_dshPert = (PfdpPert - Pfdp) / pert
    %       dPfdp_dshPert = (PfdpPert - Pfdp) / pert
    
    %%Identify power control elements (Converters and transformers included)
    iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %FUBM- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)

    %%Save the Power to control through the branch
    Pfset = branch(:,PF)/baseMVA; %Power constraint for the branch element in p.u.
    
    %%sparse matrices
    [dPfdp_dV1, dPfdp_dV2, dPfdp_dPfsh, dPfdp_dQfma, dPfdp_dBeqz,...
        dPfdp_dBeqv, dPfdp_dVtma, dPfdp_dQtma, dPfdp_dPfdp] = dPfdp_dx(branch, Yf, Yt, V, 3, 0);

    dPfdp_dV1_sp   = full(dPfdp_dV1);
    dPfdp_dV2_sp   = full(dPfdp_dV2);
    dPfdp_dPfsh_sp = full(dPfdp_dPfsh);
    dPfdp_dQfma_sp = full(dPfdp_dQfma);
    dPfdp_dBeqz_sp = full(dPfdp_dBeqz);   
    dPfdp_dBeqv_sp = full(dPfdp_dBeqv);  
    dPfdp_dVtma_sp = full(dPfdp_dVtma);
    dPfdp_dQtma_sp = full(dPfdp_dQtma);
    dPfdp_dPfdp_sp = full(dPfdp_dPfdp);
    
    %%compute numerically to compare (Finite Differences Method)
    %Voltages "from" and "to" Perturbed
    V1pf = V1p(f,:);
    V2pf = V2p(f,:);
    
    %Obtain Sf(x) and Sf(x+h) for Voltages
    Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * VV); %Sf(x)
    S1pf = V1pf .* conj(Yf * V1p); %Sf(x+h) for Va
    S2pf = V2pf .* conj(Yf * V2p); %Sf(x+h) for Vm
    
    %Obtain Pfdp(x) and Pfdp(x+h) for Voltages
    Pfdp2   = -real(Sf2)  + Pfset(:) + Kdp(:).*(Vm(f(:))  - Vmfset(:));%Pfdp(x)         | Droop Pf(x)   - Pfset = Kdp*(Vmf   - Vmfset)
    P1fdppt = -real(S1pf) + Pfset(:) + Kdp(:).*(Vm(f(:))  - Vmfset(:));%Pfdp(x+h) in Va | Droop Pf(x+h) - Pfset = Kdp*(Vmf   - Vmfset)
    P2fdppt = -real(S2pf) + Pfset(:) + Kdp(:).*(abs(V2pf) - Vmfset(:));%Pfdp(x+h) in Vm | Droop Pf(x+h) - Pfset = Kdp*(Vmf+h - Vmfset)
    
    %Voltages dPfdp_dV Finite Diff. 
    %f'(x) ~~ ( f(x+h) - f(x) ) / h 
    num_dPfdp_dV1 = full( (P1fdppt - Pfdp2) / pert );
    num_dPfdp_dV2 = full( (P2fdppt - Pfdp2) / pert );

    %FUBM extra variables Finite Diff.
    [num_dPfdp_dPfsh, num_dPfdp_dQfma, num_dPfdp_dBeqz,...
        num_dPfdp_dBeqv, num_dPfdp_dVtma, num_dPfdp_dQtma, num_dPfdp_dPfdp] = dPfdp_dxPert(baseMVA, bus, branch, Yf, Yt, V, pert, 0); %Size of each derivatives [nb, nXXxx]
    
    num_dPfdp_dPfsh_sp = full(num_dPfdp_dPfsh);
    num_dPfdp_dQfma_sp = full(num_dPfdp_dQfma);
    num_dPfdp_dBeqz_sp = full(num_dPfdp_dBeqz);
    num_dPfdp_dBeqv_sp = full(num_dPfdp_dBeqv);
    num_dPfdp_dVtma_sp = full(num_dPfdp_dVtma);
    num_dPfdp_dQtma_sp = full(num_dPfdp_dQtma);
    num_dPfdp_dPfdp_sp = full(num_dPfdp_dPfdp);
    
    %Results Comparison
    t_is(dPfdp_dV1_sp, num_dPfdp_dV1, 5, sprintf('%s - dPfdp_dV%s (sparse)', coord, vv{1}));
    t_is(dPfdp_dV2_sp, num_dPfdp_dV2, 5, sprintf('%s - dPfdp_dV%s (sparse)', coord, vv{2}));
    t_is(dPfdp_dPfsh_sp, num_dPfdp_dPfsh_sp, 5, sprintf('%s - dPfdp_dPfsh (sparse)', coord));
    t_is(dPfdp_dQfma_sp, num_dPfdp_dQfma_sp, 5, sprintf('%s - dPfdp_dQfma (sparse)', coord));
    t_is(dPfdp_dBeqz_sp, num_dPfdp_dBeqz_sp, 5, sprintf('%s - dPfdp_dBeqz (sparse)', coord));
    t_is(dPfdp_dBeqv_sp, num_dPfdp_dBeqv_sp, 5, sprintf('%s - dPfdp_dBeqv (sparse)', coord));
    t_is(dPfdp_dVtma_sp, num_dPfdp_dVtma_sp, 5, sprintf('%s - dPfdp_dVtma (sparse)', coord));
    t_is(dPfdp_dQtma_sp, num_dPfdp_dQtma_sp, 5, sprintf('%s - dPfdp_dQtma (sparse)', coord));       
    t_is(dPfdp_dPfdp_sp, num_dPfdp_dPfdp_sp, 5, sprintf('%s - dPfdp_dPfdp (sparse)', coord));

t_end;
