function [iC, iR, omegab, rc_mode] = mask_init_RC(dynamic, R, C, fbase)
%MASK_INIT_RC Pure mask initialization for dq RC branch (no set_param).

if nargin < 4 || isempty(fbase)
    fbase = 50;
end
C_safe = C;
C_safe(C_safe == 0) = 1e-6;
R_safe = R;
R_safe(R_safe == 0) = 1e6;
iC     = 1 ./ C_safe;
iR     = 1 ./ R_safe;
omegab = 2 * pi * fbase;
rc_mode = mask_parse_mode('rlrc', dynamic);
end
