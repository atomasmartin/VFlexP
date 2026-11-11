function [block_buses, block_idx, block_toM, block_fref, Lbiginv, ...
    Lssd_pu, Lssq_pu, Lmd_pu, Lmq_pu, Lfd_pu, Lkd1_pu, Lkd2_pu, Lkq1_pu, Lkq2_pu, Lkq3_pu, ...
    K_C, K_P, T_1, T_2, c_fd, c_kd1, Lppd_pu, c_kq1, c_kq2, Lppq_pu, ...
    stator_mode, use_Tmeas, use_PSS, use_LeadLag, use_KP, T_meas_safe, T_w3_safe, T_1_safe] = ...
    mask_init_MGES(blk, Lssd_pu, Lssq_pu, Lmd_pu, Lmq_pu, Lfd_pu, ...
    Lkd1_pu, Lkd2_pu, Lkq1_pu, Lkq2_pu, Lkq3_pu, ...
    K_C, K_P, T_1, T_2, dynamic_stator, T_meas, K_S1, T_w3)
%MASK_INIT_MGES Pure mask initialization for MGES (GENROU/GENSAL + ST4B + IEEEG1 + PSS2A).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'MGES');
n = max(1, numel(block_buses));

Lssd_pu = expand_n(Lssd_pu, n);
Lssq_pu = expand_n(Lssq_pu, n);
Lmd_pu  = expand_n(Lmd_pu,  n);
Lmq_pu  = expand_n(Lmq_pu,  n);
Lfd_pu  = expand_n(Lfd_pu,  n);
Lkd1_pu = expand_n(Lkd1_pu, n);
Lkd2_pu = expand_n(Lkd2_pu, n);
Lkq1_pu = expand_n(Lkq1_pu, n);
Lkq2_pu = expand_n(Lkq2_pu, n);
Lkq3_pu = expand_n(Lkq3_pu, n);
K_C     = expand_n(K_C,     n);
K_P     = expand_n(K_P,     n);
T_1     = expand_n(T_1,     n);
T_2     = expand_n(T_2,     n);

Lbiginv = zeros(9*n, 9*n);
for i = 1:n
    j = 9*(i-1) + 1;
    Ldq0  = [Lssd_pu(i), 0, 0; 0, Lssq_pu(i), 0; 0, 0, 0.0001];
    Ladq0 = [Lmd_pu(i), Lmd_pu(i), Lmd_pu(i), 0, 0, 0; ...
             0, 0, 0, Lmq_pu(i), Lmq_pu(i), Lmq_pu(i); ...
             zeros(1, 6)];
    Lk = [Lmd_pu(i)*ones(3,3), zeros(3,3); zeros(3,3), Lmq_pu(i)*ones(3,3)] + ...
         diag([Lfd_pu(i), Lkd1_pu(i), Lkd2_pu(i), Lkq1_pu(i), Lkq2_pu(i), Lkq3_pu(i)]);
    L = [-Ldq0, Ladq0; -Ladq0', Lk];
    Lbiginv(j:(j+8), j:(j+8)) = inv(L);
end

% Subtransient coefficients for order 4 (RMS4) algebraic stator
denom_d = Lfd_pu .* Lkd1_pu + Lmd_pu .* (Lfd_pu + Lkd1_pu);
denom_d(denom_d == 0) = 1;
c_fd    = Lmd_pu .* Lkd1_pu ./ denom_d;
c_kd1   = Lmd_pu .* Lfd_pu  ./ denom_d;
Lppd_pu = Lssd_pu - Lmd_pu + Lmd_pu .* Lfd_pu .* Lkd1_pu ./ denom_d;

denom_q = Lkq1_pu .* Lkq2_pu + Lmq_pu .* (Lkq1_pu + Lkq2_pu);
denom_q(denom_q == 0) = 1;
c_kq1   = Lmq_pu .* Lkq2_pu ./ denom_q;
c_kq2   = Lmq_pu .* Lkq1_pu ./ denom_q;
Lppq_pu = Lssq_pu - Lmq_pu + Lmq_pu .* Lkq1_pu .* Lkq2_pu ./ denom_q;

dynStatVal  = mask_read_prm_override(blk, 'MGES', 'dynamic_stator', dynamic_stator);
stator_mode = mask_parse_mode('stator', dynStatVal);

use_Tmeas   = double(~all(T_meas == 0));
use_PSS     = double(~(all(K_S1 == 0) || all(T_w3 == 0)));
use_LeadLag = double(~all(T_1 == 0));
use_KP      = double(~all(K_P == 0));

T_meas_safe = T_meas;
T_meas_safe(T_meas_safe == 0) = 0.01;
T_w3_safe = T_w3;
T_w3_safe(T_w3_safe == 0) = 1;
T_1_safe = T_1;
T_1_safe(T_1_safe == 0) = 1;
end

function v = expand_n(x, n)
if isscalar(x)
    v = x * ones(1, n);
else
    v = reshape(x, 1, []);
end
end
