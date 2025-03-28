function fun_initialise_0_callback(gcb)
    bl_int = find_system([gcb '/Dynamic'], 'LookUnderMasks', 'on', 'BlockType', 'Integrator');
    for i=1:length(bl_int)
        set_param(bl_int{i}, 'InitialCondition', '0')
    end
end