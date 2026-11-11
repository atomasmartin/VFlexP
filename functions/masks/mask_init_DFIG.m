function [block_buses, block_idx, block_toM, block_fref, Lbiginv, ...
    Ls_pu, Lr_pu, Lm_pu, La_pu, Lss_pu, Lrr_pu, sigma, omega_base] = ...
    mask_init_DFIG(blk, FBASE, Ls_pu, Lr_pu, Lm_pu, La_pu)
%MASK_INIT_DFIG Pure mask initialization for DFIG.

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'DFIG');
n = max(1, numel(block_buses));

if isscalar(Ls_pu), Ls_pu = Ls_pu * ones(1, n); end
if isscalar(Lr_pu), Lr_pu = Lr_pu * ones(1, n); end
if isscalar(Lm_pu), Lm_pu = Lm_pu * ones(1, n); end
if isscalar(La_pu), La_pu = La_pu * ones(1, n); end

Lss_pu = Ls_pu + Lm_pu;
Lrr_pu = Lr_pu + Lm_pu;

Lbiginv = zeros(6*n, 6*n);
for i = 1:n
    j = 6*(i-1) + 1;
    L = [Lss_pu(i), 0, Lm_pu(i), 0, 0, 0; ...
         0, Lss_pu(i), 0, Lm_pu(i), 0, 0; ...
         Lm_pu(i), 0, Lrr_pu(i), 0, 0, 0; ...
         0, Lm_pu(i), 0, Lrr_pu(i), 0, 0; ...
         0, 0, 0, 0, La_pu(i), 0; ...
         0, 0, 0, 0, 0, La_pu(i)];
    Lbiginv(j:(j+5), j:(j+5)) = inv(L);
end

sigma      = 1 - Lm_pu.^2 ./ (Lss_pu .* Lrr_pu);
omega_base = FBASE * 2 * pi;
end
