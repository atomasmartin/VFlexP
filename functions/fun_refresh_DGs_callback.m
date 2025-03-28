function fun_refresh_DGs_callback(gcb)
% AC grid
model = bdroot; % current model
block = gcb;
bl = [find_system([block '/Dynamic/Generators'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on'); ...
    find_system([block '/Dynamic/Generators_DC'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on')];
bl_names = erase(bl,[block '/Dynamic/Generators/']);
bl_names = erase(bl_names,[block '/Dynamic/Generators_DC/']);

if ~bdIsLibrary(model)
    maskObj = Simulink.Mask.get(gcb);
    tableControl = maskObj.getDialogControl('gen_data');
    options = {tableControl.Columns(:).TypeOptions};
    idx = find(~cellfun(@isempty, options));
    tableControl.Columns(1, idx).TypeOptions = bl_names;
end
% Interfacing converters
bl = find_system([block '/Dynamic/DCAC converters'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on');
bl_names = erase(bl,[block '/Dynamic/DCAC converters/']);
bl_names = erase(bl_names,[block '/Dynamic/DCAC converters']);

if ~bdIsLibrary(model)
    maskObj = Simulink.Mask.get(gcb);
    tableControl = maskObj.getDialogControl('intconv_data');
    options = {tableControl.Columns(:).TypeOptions};
    idx = find(~cellfun(@isempty, options));
    tableControl.Columns(1, idx).TypeOptions = bl_names;
end
end