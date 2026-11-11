function [block_buses, block_idx, block_toM, block_fref] = mask_init_IG(blk)
%MASK_INIT_IG Pure mask initialization for IG (infinite generator / ideal voltage source).

[block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, 'IG');
end
