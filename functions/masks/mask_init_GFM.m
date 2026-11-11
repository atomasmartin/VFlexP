function [block_buses, block_idx, block_toM, block_fref, vcc_mode] = mask_init_GFM(blk, dynamic_VCC)
%MASK_INIT_GFM Pure mask initialization for GFM (droop grid-forming converter).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'GFM');
vccVal   = mask_read_prm_override(blk, 'GFM', 'dynamic_VCC', dynamic_VCC);
vcc_mode = mask_parse_mode('vcc', vccVal);
end
