function [invL, omegab, rl_mode] = mask_init_RL(dynamic, L, fbase)
%MASK_INIT_RL Pure mask initialization for dq RL branch (no set_param).

if nargin < 3 || isempty(fbase)
    fbase = 50;
end
L_safe = L;
L_safe(L_safe == 0) = 1e-6;
invL   = 1 ./ L_safe;
omegab = 2 * pi * fbase;
rl_mode = mask_parse_mode('rlrc', dynamic);
end
