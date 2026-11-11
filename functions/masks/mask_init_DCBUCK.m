function [block_buses, block_Pg, block_toM, block_idx] = mask_init_DCBUCK(blk)
%MASK_INIT_DCBUCK Pure mask initialization for DCBUCK DC generator.

[block_buses, block_Pg, block_toM, block_idx] = mask_dc_unit_routing(blk, 'DCBUCK');
end
