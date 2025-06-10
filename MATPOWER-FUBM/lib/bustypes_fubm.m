function [ref, pv, pq] = bustypes_fubm(bus, gen, branch)
%BUSTYPES_FUBM   Builds index lists for each type of bus (REF, PV, PQ).
%   [REF, PV, PQ] = BUSTYPES(BUS, GEN, BRANCH)
%   Generators with "out-of-service" status are treated as PQ buses with
%   zero generation (regardless of Pg/Qg values in gen). Expects BUS and
%   GEN have been converted to use internal consecutive bus numbering.

%   Voltage controlled buses are considered as PV buses.

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM- Original: idx_brch
%% get generator status
% bus_gen_status = zeros(size(bus, 1), 1);
% bus_gen_status(gen(:, GEN_BUS)) = gen(:, GEN_STATUS) > 0;
nb = size(bus, 1);
ng = size(gen, 1);
Cg = sparse(gen(:, GEN_BUS), (1:ng)', gen(:, GEN_STATUS) > 0, nb, ng);  %% gen connection matrix

                                        %% element i, j is 1 if, generator j at bus i is ON
bus_gen_status = Cg * ones(ng, 1);      %% number of generators at each bus that are ON
%% Check if is FUBM
%%AAB----------------------------------------------------------------------
if size(branch,2) < ALPH3
    fubm = 0;
else
    fubm = 1;
end
%%-------------------------------------------------------------------------
if fubm
    %% Check if there are elements that control the Voltage magnitude
    %%AAB----------------------------------------------------------------------
    %find elements that Control Vm and their buses
    iVt_ctrl = find((branch(:,BR_STATUS)~=0) & (branch(:,VT_SET)~=0)); %location of the elements that control Vt
    iVf_ctrl = find((branch(:,BR_STATUS)~=0) & (branch(:,VF_SET)~=0) & ( branch(:, CONV)== 2 | (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) )); %location of the elements that control Vf
    
    Vt_ctrl_bus = branch(iVt_ctrl, T_BUS); %Locate the  "to"  bus of the elements with voltage  "to"  control 
    Vf_ctrl_bus = branch(iVf_ctrl, F_BUS); %Locate the "from" bus of the elements with voltage "from" control
    %%-------------------------------------------------------------------------

    %% Change type of bus for Voltage controlled elements
    %%AAB----------------------------------------------------------------------
    bus(Vf_ctrl_bus, BUS_TYPE) = 2; %Change the PQ bus to a PV bus. 
    bus(Vt_ctrl_bus, BUS_TYPE) = 2; %Change the PQ bus to a PV bus.
    %%-------------------------------------------------------------------------

    %% Correct the bus type when tap Changer Voltage controlled is included
    %%AAB---------------------------------------------------------------------
    %%Get branch status
    nl = size(branch, 1); %number of lines
    vctrl_t=zeros(nl,1);  %vector with the same size as existing lines for Voltage control "to"
    vctrl_f=zeros(nl,1);  %vector with the same size as existing lines for Voltage control "from"
    vctrl_t(find(branch(:,VT_SET))) = 1; %Sets a value of 1 when the line controls voltage "to" side.
    vctrl_f(find(branch(:,VF_SET))) = 1; %Sets a value of 1 when the line controls voltage "from" side.
    Cvclt = sparse(branch(:, T_BUS), (1:nl)', vctrl_t > 0, nb, nl);  %% voltage control branch connection matrix
    Cvclf = sparse(branch(:, F_BUS), (1:nl)', vctrl_f > 0, nb, nl);  %% voltage control branch connection matrix

                                                    %%element i, j is 1 if, branch j at bus i is ON
    bus_branch_status_t = Cvclt * ones(nl, 1);      %% number of voltage control branches at each bus that are ON
    bus_branch_status_f = Cvclf * ones(nl, 1);      %% number of voltage control branches at each bus that are ON
    %%-------------------------------------------------------------------------
end
%% form index lists for slack, PV, and PQ buses (Based On Generation)
%%AAB----------------------------------------------------------------------
if fubm
    ref = find(bus(:, BUS_TYPE) == REF & (bus_gen_status | bus_branch_status_t | bus_branch_status_f));   %% reference bus index 
    pv  = find(bus(:, BUS_TYPE) == PV  & (bus_gen_status | bus_branch_status_t | bus_branch_status_f));   %% PV bus indices
    pq  = find(bus(:, BUS_TYPE) == PQ | ~(bus_gen_status | bus_branch_status_t | bus_branch_status_f));   %% PQ bus indices
else
    ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   %% reference bus index
    pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   %% PV bus indices
    pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);    
end    
%%-------------------------------------------------------------------------
%% pick a new reference bus if for some reason there is none (may have been shut down)
if isempty(ref)
%     if isempty(pv)
%         %% no PV bus left to convert to reference bus
%     else
        ref = pv(1);    %% use the first PV bus
        pv(1) = [];     %% delete it from PV list
%     end
end
