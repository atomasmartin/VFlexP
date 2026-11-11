function [block_buses, block_idx, block_toM, block_fref, Lbiginv, ...
    Lssd_pu, Lssq_pu, Lmd_pu, Lmq_pu, Lfd_pu, Lkd1_pu, Lkd2_pu, Lkq1_pu, Lkq2_pu, Lkq3_pu] = ...
    mask_init_SC(blk, Lssd_pu, Lssq_pu, Lmd_pu, Lmq_pu, Lfd_pu, Lkd1_pu, Lkd2_pu, Lkq1_pu, Lkq2_pu, Lkq3_pu)
%MASK_INIT_SC Pure mask initialization for SC (synchronous condenser).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'SC');
n = max(1, numel(block_buses));

if isscalar(Lssd_pu), Lssd_pu = Lssd_pu * ones(1, n); end
if isscalar(Lssq_pu), Lssq_pu = Lssq_pu * ones(1, n); end
if isscalar(Lmd_pu),  Lmd_pu  = Lmd_pu  * ones(1, n); end
if isscalar(Lmq_pu),  Lmq_pu  = Lmq_pu  * ones(1, n); end
if isscalar(Lfd_pu),  Lfd_pu  = Lfd_pu  * ones(1, n); end
if isscalar(Lkd1_pu), Lkd1_pu = Lkd1_pu * ones(1, n); end
if isscalar(Lkd2_pu), Lkd2_pu = Lkd2_pu * ones(1, n); end
if isscalar(Lkq1_pu), Lkq1_pu = Lkq1_pu * ones(1, n); end
if isscalar(Lkq2_pu), Lkq2_pu = Lkq2_pu * ones(1, n); end
if isscalar(Lkq3_pu), Lkq3_pu = Lkq3_pu * ones(1, n); end

Lbiginv = zeros(9*n, 9*n);
for i = 1:n
    j = 9*(i-1) + 1;
    Ldq0  = [Lssd_pu(i), 0, 0; 0, Lssq_pu(i), 0; 0, 0, 0.1];
    Ladq0 = [Lmd_pu(i), Lmd_pu(i), Lmd_pu(i), 0, 0, 0; ...
             0, 0, 0, Lmq_pu(i), Lmq_pu(i), Lmq_pu(i); ...
             zeros(1, 6)];
    Lk = [Lmd_pu(i)*ones(3,3), zeros(3,3); zeros(3,3), Lmq_pu(i)*ones(3,3)] + ...
         diag([Lfd_pu(i), Lkd1_pu(i), Lkd2_pu(i), Lkq1_pu(i), Lkq2_pu(i), Lkq3_pu(i)]);
    L = [-Ldq0, Ladq0; -Ladq0', Lk];
    Lbiginv(j:(j+8), j:(j+8)) = inv(L);
end
end
