function fun_load_complexity_grid(gcb)
    [Tbus, ~, Tline, TbusDC, TlineDC, ~, ~, ~, ~, TgenDC] = fun_read_mask_tables(gcb);
    fun_set_complexity_blocks_grid(gcb, Tline, Tbus, TbusDC, TlineDC, TgenDC)
end