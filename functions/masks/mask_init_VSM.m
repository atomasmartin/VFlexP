function [block_buses, block_idx, block_toM, block_fref, vcc_mode] = mask_init_VSM(blk, dynamic_VCC)
%MASK_INIT_VSM Pure mask initialization for VSM (virtual synchronous machine).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'VSM');
vccVal   = mask_read_prm_override(blk, 'VSM', 'dynamic_VCC', dynamic_VCC);
vcc_mode = mask_parse_mode('vcc', vccVal);
end
