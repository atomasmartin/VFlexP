function fun_fill_tables_in_mask(gcb, Tbus, Tgen, Tline, Tintconv, Tparam)
try
    block = gcb;
    % Bus table
    str = '{';
    Rows = height(Tbus);
    Columns = width(Tbus);
    for i = 1:Rows
        for ii = 1:Columns
            if ismember(ii, [8 9])
                if iscell(Tbus{i,ii})
                    str = [str char(" '") char(num2str(Tbus{i,ii}{1}, '%.10f')) char("' ")];
                else
                    str = [str char(" '") char(num2str(Tbus{i,ii}, '%.10f')) char("' ")];
                end
            else
                str = [str char(" '") char(string(Tbus{i,ii})) char("' ")];
            end
        end
        if i < Rows
            str = [str,';'];
        end
    end
    str = [str,'}'];
    set_param(block, "bus_data", str)
    % Gen table
    str = '{';
    Rows = height(Tgen);
    Columns = width(Tgen);
    for i = 1:Rows
        for ii = 1:Columns
            if ismember(ii, [3 4])
                if iscell(Tgen{i,ii})
                    str = [str char(" '") char(num2str(Tgen{i,ii}{1}, '%.10f')) char("' ")];
                else
                    str = [str char(" '") char(num2str(Tgen{i,ii}, '%.10f')) char("' ")];
                end
            else
                str = [str char(" '") char(string(Tgen{i,ii})) char("' ")];
            end
        end
        if i < Rows
            str = [str,';'];
        end
    end
    str = [str,'}'];
    set_param(block, "gen_data", str)
    % Line table
    str = '{';
    Rows = height(Tline);
    Columns = width(Tline);
    for i = 1:Rows
        for ii = 1:Columns
            str = [str char(" '") char(string(Tline{i,ii})) char("' ")];
        end
        if i < Rows
            str = [str,';'];
        end
    end
    str = [str,'}'];
    set_param(block, "line_data", str)
    % Interface converters table
    if ~isempty(Tintconv)
        str = '{';
        Rows = height(Tintconv);
        Columns = width(Tintconv);
        for i = 1:Rows
            for ii = 1:Columns
                str = [str char(" '") char(string(Tintconv{i,ii})) char("' ")];
            end
            if i < Rows
                str = [str,';'];
            end
        end
        str = [str,'}'];
        set_param(block, "intconv_data", str)
    end

    if nargin > 5
        % Parameter variation table
        if ~isempty(Tparam)
            str = '{';
            Rows = height(Tparam);
            Columns = width(Tparam);
            for i = 1:Rows
                for ii = 1:Columns
                    str = [str char(" '") char(string(Tparam{i,ii})) char("' ")];
                end
                if i < Rows
                    str = [str,';'];
                end
            end
            str = [str,'}'];
            set_param(block, "param_change_data", str)
        end
    end
    fun_apply_changes_block(gcb)
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end