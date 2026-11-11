function [block_buses, block_idx, block_toM, block_fref, iTf, use_Tf, cc_mode, pll_mode] = ...
    mask_init_GFL(blk, Tf, dynamic_CC, dynamic_PLL)
%MASK_INIT_GFL Pure mask initialization for GFL (grid-following converter).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'GFL');
[iTf, use_Tf] = mask_init_opt_filter(Tf);

ccVal    = mask_read_prm_override(blk, 'GFL', 'dynamic_CC',  dynamic_CC);
pllVal   = mask_read_prm_override(blk, 'GFL', 'dynamic_PLL', dynamic_PLL);
cc_mode  = mask_parse_mode('onoff', ccVal);
pll_mode = mask_parse_mode('onoff', pllVal);
end
