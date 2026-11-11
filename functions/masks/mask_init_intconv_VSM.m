function [block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx, ...
    vcc_mode, use_wof_dc, use_kiw] = mask_init_intconv_VSM(blk, dynamic_VCC, Kd_AC, Ki_w)
%MASK_INIT_INTCONV_VSM Pure mask initialization for intconv_VSM DC/AC converter.

[block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx] = ...
    mask_ic_routing(blk, 'intconv_VSM');
vccVal     = mask_read_prm_override(blk, 'intconv_VSM', 'dynamic_VCC', dynamic_VCC);
vcc_mode   = mask_parse_mode('vcc', vccVal);
use_wof_dc = double(~all(Kd_AC == 0));
use_kiw    = double(~all(Ki_w == 0));
end
