function fun_save_table_callback(gcb, excel_file)
try
    block = gcb;
    maskObj = Simulink.Mask.get(block);

    % Function to safely prepare table from raw data
    % defined at the end of the file.

    % Bus data
    var_names1 = string({maskObj.getDialogControl('bus_data').Columns(:).Name});
    raw_bus_data = eval(get_param(block, 'bus_data'));
    table_values1 = local_prepare_table(raw_bus_data, var_names1);

    % Line data
    var_names2 = string({maskObj.getDialogControl('line_data').Columns(:).Name});
    raw_line_data = eval(get_param(block, 'line_data'));
    if ~isempty(raw_line_data)
        table_values2 = local_prepare_table(raw_line_data, var_names2);
    end

    % Gen data
    var_names3 = string({maskObj.getDialogControl('gen_data').Columns(:).Name});
    raw_gen_data = eval(get_param(block, 'gen_data'));
    table_values3 = local_prepare_table(raw_gen_data, var_names3);

    % IntConv data
    var_names4 = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
    raw_intconv = eval(get_param(block, 'intconv_data'));

    if nargin == 1
        [table_file, path] = uiputfile('.xlsx');
        if isequal(table_file, 0) || isequal(path, 0)
            return;
        end
    else
        [path, table_file, ~] = fileparts(excel_file);
    end

    save_path = fullfile(path, table_file);

    writetable(table_values1, save_path, 'Sheet', "Bus data")
    if ~isempty(raw_line_data)
        writetable(table_values2, save_path, 'Sheet', "Line data")
    end
    writetable(table_values3, save_path, 'Sheet', "Gen data")

    if ~isempty(raw_intconv)
        table_values4 = local_prepare_table(raw_intconv, var_names4);
        writetable(table_values4, save_path, 'Sheet', "IntConv data")
    end

    % Fb_data (System parameters)
    raw_fb_data = eval(get_param(block, 'Fb_data'));
    if ~isempty(raw_fb_data)
        var_names5 = string({maskObj.getDialogControl('Fb_data').Columns(:).Name});
        table_values5 = local_prepare_table(raw_fb_data, var_names5);

        % Add Sb column based on last column
        table_values5.Sb = table_values5{:, end};
        sb_val = evalin('base', get_param(gcb, 'Sb'));

        % Handle Sb(1) update (case for numeric or cell column)
        if iscell(table_values5.Sb)
            table_values5.Sb{1} = sb_val;
        else
            table_values5.Sb(1) = sb_val;
        end

        movevars(table_values5, 'Sb', 'Before', 'Area')
        writetable(table_values5, save_path, 'Sheet', "System parameters")
    end

    % Parameter change data
    raw_param_data = eval(get_param(block, 'param_change_data'));
    if ~isempty(raw_param_data)
        var_names6 = string({maskObj.getDialogControl('param_change_data').Columns(:).Name});
        table_values6 = local_prepare_table(raw_param_data, var_names6);
        writetable(table_values6, save_path, 'Sheet', "Parameter variation")
    end

catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end

function table_out = local_prepare_table(data, var_names)
% Helper to convert raw mask data (cell or array) to a table with numeric types where possible
if isempty(data)
    table_out = table.empty(0, length(var_names));
    table_out.Properties.VariableNames = var_names;
    return;
end

if ~iscell(data)
    if isnumeric(data)
        table_out = array2table(data, 'VariableNames', var_names);
        return;
    else
        data = cellstr(data);
    end
end

% Try to convert numeric-like strings to double for each element
for i = 1:numel(data)
    if ischar(data{i}) || isstring(data{i})
        num = str2double(data{i});
        if ~isnan(num)
            data{i} = num;
        end
    end
end

table_out = cell2table(data, 'VariableNames', var_names);

% Post-process: simplify cell columns to numeric arrays if they are purely numeric
for c = 1:width(table_out)
    col_data = table_out{:, c};
    if iscell(col_data) && all(cellfun(@(x) isnumeric(x) && isscalar(x), col_data))
        try
            table_out.(table_out.Properties.VariableNames{c}) = cell2mat(col_data);
        catch
        end
    end
end
end
