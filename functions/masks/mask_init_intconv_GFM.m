function [block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx, vcc_mode] = ...
    mask_init_intconv_GFM(blk, dynamic_VCC)
%MASK_INIT_INTCONV_GFM Pure mask initialization for intconv_GFM DC/AC converter.

[block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx] = ...
    mask_ic_routing(blk, 'intconv_GFM');
vccVal   = mask_read_prm_override(blk, 'intconv_GFM', 'dynamic_VCC', dynamic_VCC);
vcc_mode = mask_parse_mode('vcc', vccVal);
end
