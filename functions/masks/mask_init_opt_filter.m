function [iTf, use_Tf] = mask_init_opt_filter(Tf)
%MASK_INIT_OPT_FILTER Pure mask initialization for optional LPF / WOF sub-blocks (PC, VI).

use_Tf = double(~all(Tf == 0));
Tf_safe = Tf;
Tf_safe(Tf_safe == 0) = 1;
iTf = 1 ./ Tf_safe;
end
