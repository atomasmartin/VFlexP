function [block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx, ...
    iTf, use_Tf, cc_mode, pll_mode] = mask_init_intconv_GFL(blk, Tf, dynamic_CC, dynamic_PLL)
%MASK_INIT_INTCONV_GFL Pure mask initialization for intconv_GFL DC/AC converter.

[block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx] = ...
    mask_ic_routing(blk, 'intconv_GFL');
[iTf, use_Tf] = mask_init_opt_filter(Tf);

ccVal    = mask_read_prm_override(blk, 'intconv_GFL', 'dynamic_CC',  dynamic_CC);
pllVal   = mask_read_prm_override(blk, 'intconv_GFL', 'dynamic_PLL', dynamic_PLL);
cc_mode  = mask_parse_mode('onoff', ccVal);
pll_mode = mask_parse_mode('onoff', pllVal);
end
