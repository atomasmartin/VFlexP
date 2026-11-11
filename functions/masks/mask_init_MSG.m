function [block_buses, block_idx, block_toM, block_fref, Lbiginv, ...
    Lssd_pu, Lssq_pu, Lmd_pu, Lfd_pu, La_pu, omega_base] = ...
    mask_init_MSG(blk, FBASE, Lssd_pu, Lssq_pu, Lmd_pu, Lfd_pu, La_pu)
%MASK_INIT_MSG Pure mask initialization for MSG (motor-synchronous generator).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'MSG');
n = max(1, numel(block_buses));

if isscalar(Lssd_pu), Lssd_pu = Lssd_pu * ones(1, n); end
if isscalar(Lssq_pu), Lssq_pu = Lssq_pu * ones(1, n); end
if isscalar(Lmd_pu),  Lmd_pu  = Lmd_pu  * ones(1, n); end
if isscalar(Lfd_pu),  Lfd_pu  = Lfd_pu  * ones(1, n); end
if isscalar(La_pu),   La_pu   = La_pu   * ones(1, n); end

Lbiginv = zeros(6*n, 6*n);
for i = 1:n
    j = 6*(i-1) + 1;
    Ldq0  = [Lssd_pu(i), 0, 0; 0, Lssq_pu(i), 0; 0, 0, Lssd_pu(i) - Lmd_pu(i)];
    Ladq0 = [Lmd_pu(i), 0, 0; 0, 0, 0; zeros(1, 3)];
    Lk    = [Lmd_pu(i), 0, 0; zeros(1, 3); zeros(1, 3)] + diag([Lfd_pu(i), La_pu(i), La_pu(i)]);
    L     = [Ldq0, Ladq0; Ladq0', Lk];
    Lbiginv(j:(j+5), j:(j+5)) = inv(L);
end

omega_base = FBASE * 2 * pi;
end
