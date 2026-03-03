function fun_refresh_DGs_callback(gcb)
% AC grid
model = bdroot; % current model
block = gcb;
bl = [find_system([block '/AC Generators'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on'); ...
    find_system([block '/DC Generators'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on')];
bl_names = erase(bl,[block '/AC Generators/']);
bl_names = erase(bl_names,[block '/DC Generators/']);

if ~bdIsLibrary(model)
    maskObj = Simulink.Mask.get(gcb);
    tableControl = maskObj.getDialogControl('gen_data');
    options = {tableControl.Columns(:).TypeOptions};
    idx = find(~cellfun(@isempty, options));
    tableControl.Columns(1, idx).TypeOptions = bl_names;
end
% Interfacing converters
bl = find_system([block '/DCAC Converters'], ...
    'IncludeCommented', 'on', 'SearchDepth', '1', ...
    'LookUnderMasks', 'on', 'Mask','on');
bl_names = erase(bl,[block '/DCAC Converters/']);
bl_names = erase(bl_names,[block '/DCAC Converters']);

if ~bdIsLibrary(model)
    maskObj = Simulink.Mask.get(gcb);
    tableControl = maskObj.getDialogControl('intconv_data');
    options = {tableControl.Columns(:).TypeOptions};
    idx = find(~cellfun(@isempty, options));
    tableControl.Columns(1, idx).TypeOptions = bl_names;
end
end