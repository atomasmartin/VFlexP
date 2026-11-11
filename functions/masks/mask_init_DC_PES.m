function [block_buses, block_Pg, block_toM, block_idx, use_wof] = mask_init_DC_PES(blk, Kd_dc)
%MASK_INIT_DC_PES Pure mask initialization for DC_PES DC primary energy source.

[block_buses, block_Pg, block_toM, block_idx] = mask_dc_unit_routing(blk, 'DC_PES');
use_wof = double(~all(Kd_dc == 0));
end
