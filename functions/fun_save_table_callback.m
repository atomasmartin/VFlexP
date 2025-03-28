function fun_save_table_callback(gcb)
try
    block = gcb;
    maskObj = Simulink.Mask.get(block);
    var_names1 = string({maskObj.getDialogControl('bus_data').Columns(:).Name});
    table_values1 = array2table(string(eval(get_param(block, 'bus_data'))), ...
        "VariableNames", var_names1);
    var_names2 = string({maskObj.getDialogControl('line_data').Columns(:).Name});
    line_table_string = string(eval(get_param(block, 'line_data')));
    if ~isempty(line_table_string)
        table_values2 = array2table(line_table_string, "VariableNames", var_names2);
    end
    var_names3 = string({maskObj.getDialogControl('gen_data').Columns(:).Name});
    table_values3 = array2table(string(eval(get_param(block, 'gen_data'))), ...
        "VariableNames", var_names3);
    var_names4 = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
    str_intconv = string(eval(get_param(block, 'intconv_data')));


    [table_file, path] = uiputfile('.xlsx');
    writetable(table_values1, fullfile(path, table_file), 'Sheet', "Bus data")
    if ~isempty(line_table_string)
        writetable(table_values2, fullfile(path, table_file), 'Sheet', "Line data")
    end
    writetable(table_values3, fullfile(path, table_file), 'Sheet', "Gen data")
    if ~isempty(str_intconv)
        table_values4 = array2table(str_intconv, "VariableNames", var_names4);
        writetable(table_values4, fullfile(path, table_file), 'Sheet', "IntConv data")
    end
    if ~isempty(string(eval(get_param(block, 'Fb_data'))))
        var_names5 = string({maskObj.getDialogControl('Fb_data').Columns(:).Name});
        table_values5 = array2table(string(eval(get_param(block, 'Fb_data'))), ...
            "VariableNames", var_names5);
        table_values5.Sb = table_values5{:,end};
        table_values5.Sb(1) = evalin('base', get_param(gcb, 'Sb'));
        movevars(table_values5, 'Sb', 'Before', 'Area')
        writetable(table_values5, fullfile(path, table_file), 'Sheet', "System parameters")
    end
    if ~isempty(string(eval(get_param(block, 'param_change_data'))))
        var_names6 = string({maskObj.getDialogControl('param_change_data').Columns(:).Name});
        table_values6 = array2table(string(eval(get_param(block, 'param_change_data'))), ...
            "VariableNames", var_names6);
        writetable(table_values6, fullfile(path, table_file), 'Sheet', "Parameter variation")
    end
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end