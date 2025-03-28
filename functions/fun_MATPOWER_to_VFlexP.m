function fun_MATPOWER_to_VFlexP(mpc, gcb, write_branch)
try

    block = gcb;
    set_param(block, 'Sb', num2str(mpc.baseMVA));
    maskObj = Simulink.Mask.get(block);

    % Bus table
    var_names = string({maskObj.getDialogControl('bus_data').Columns(:).Name});
    Tbus = array2table(mpc.bus, "VariableNames", var_names(1:end-4));
    Tbus.Complexity_bus = strings(height(Tbus),1);
    Tbus.Complexity_bus(:) = 'static omega';
    Tbus.Complexity_load = strings(height(Tbus),1);
    Tbus.Complexity_load(:) = 'static omega';
    Tbus.Complexity_bus_LF = Tbus.Complexity_bus;
    Tbus.Complexity_load_LF = Tbus.Complexity_bus;

    % Gen table
    var_names = string({maskObj.getDialogControl('gen_data').Columns(:).Name});
    mpc.gen = [mpc.gen(:,1) mpc.gen];
    Tgen = array2table(mpc.gen, "VariableNames", var_names);
    Tgen.type = repmat("GFr", height(Tgen), 1); % by default, grid-forming converter with droop control

    % Line table
    if ~isempty(mpc.branch)
        var_names = string({maskObj.getDialogControl('line_data').Columns(:).Name});
        if size(mpc.branch, 2) > 13 % FUBM formulation
            idx_lines = mpc.branch(:,26) == 0;
            idx_lines_AC = mpc.branch(idx_lines,4) ~= 0;
            idx_lines_DC = mpc.branch(idx_lines,4) == 0;
            idx_branch_DC =  mpc.branch(:,4) == 0 & idx_lines;
            DC_buses = unique([mpc.branch(idx_branch_DC,1); mpc.branch(idx_branch_DC,2)]);
            idx_buses_DC = ismember(Tbus.bus_i, DC_buses);
            Tbus.Complexity_bus(idx_buses_DC) = 'DC static omega';
            Tbus.Complexity_bus_LF(idx_buses_DC) = 'DC static omega';
            Tline = array2table(mpc.branch(idx_lines,(1:13)), "VariableNames", var_names(1:end-2));
            Tline.Complexity = strings(height(Tline),1);
            Tline.Complexity_PF = strings(height(Tline),1);
            Tline.Complexity(idx_lines_DC) = 'DC static omega';
            Tline.Complexity(idx_lines_AC) = 'static omega';
            Tline.Complexity_PF(idx_lines_DC) = 'DC static omega';
            Tline.Complexity_PF(idx_lines_AC) = 'static omega';
        else
            Tline = array2table(mpc.branch, "VariableNames", var_names(1:end-2));
            Tline.Complexity = strings(height(Tline),1);
            Tline.Complexity(:) = 'static omega';
            Tline.Complexity_LF = Tline.Complexity;
        end

        % Interface converters table
        var_names = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
        if size(mpc.branch, 2) > 13 % FUBM formulation
            idx_intconv = mpc.branch(:,26) ~= 0;
            Nintconv = sum(idx_intconv);
            if Nintconv > 0
                type = strings(Nintconv, 1);
                type(:) = "DC_GFl";
                Tintconv = [table(type) array2table(mpc.branch(idx_intconv,:), "VariableNames", var_names(2:end))];
            else
                Tintconv = array2table([], "VariableNames", var_names(1:end));
            end
        else
            Tintconv = [];
        end
    else
        Tintconv = []; Tline = [];
    end

    % If there are variables in the mask, do not overwrite them
    Sgen0 = string(eval(get_param(gcb, 'gen_data')));
    Sline0 = string(eval(get_param(gcb, 'line_data')));
    Sbus0 = string(eval(get_param(gcb, 'bus_data')));
    Sintconv0 = string(eval(get_param(gcb, 'intconv_data')));

    % Transform to num, the NaN elements were strings, do not replace them
    Agen0 = str2double(Sgen0);
    Aline0 = str2double(Sline0);
    Abus0 = str2double(Sbus0);
    Aintconv0 = str2double(Sintconv0);

    if height(Tgen) == size(Sgen0, 1)
        for i=1:width(Tgen)
            nan_elements = isnan(Agen0(:,i));
            if any(nan_elements)
                Tgen = convertvars(Tgen,i,'string'); % Transform column to strings
                Tgen{nan_elements,i} = Sgen0(nan_elements,i); % Keep it
            end
        end
    end

    if height(Tline) == (size(Sline0, 1) + size(Sintconv0, 1))
        for i=1:width(Tline)
            nan_elements = isnan(Aline0(:,i));
            if any(nan_elements)
                Tline = convertvars(Tline,i,'string'); % Transform column to strings
                Tline{nan_elements,i} = Sline0(nan_elements,i); % Keep it
            end
        end
    end

    if height(Tbus) == size(Sbus0, 1)
        for i=1:width(Tbus)
            nan_elements = isnan(Abus0(:,i));
            if any(nan_elements)
                Tbus = convertvars(Tbus,i,'string'); % Transform column to strings
                Tbus{nan_elements,i} = Sbus0(nan_elements,i); % Keep it
            end
        end
    end

    if height(Tintconv) == size(Sintconv0, 1)
        for i=1:width(Tintconv)
            nan_elements = isnan(Aintconv0(:,i));
            if any(nan_elements)
                Tintconv = convertvars(Tintconv,i,'string'); % Transform column to strings
                Tintconv{nan_elements,i} = Sintconv0(nan_elements,i); % Keep it
            end
        end
    end

    switch write_branch
        case 1
            % nothing to keep
        case 2
            [~, ~, ~, ~, ~, Tintconv, ~, Tline] = fun_read_mask_tables(gcb);
            % Interface converters table
            var_names = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
            intconv = mpc.branch_complete(mpc.branch_complete(:,26) ~= 0, :);
            Tintconv(:,2:end) = array2table(intconv, "VariableNames", var_names(2:end));
        case 3
            var_names = string({maskObj.getDialogControl('line_data').Columns(:).Name});
            lines = eval(get_param(block, 'line_data'));
            Tline = array2table(lines, "VariableNames", var_names);
        case 4
            var_names = string({maskObj.getDialogControl('line_data').Columns(:).Name});
            lines = eval(get_param(block, 'line_data'));
            Tline = array2table(lines, "VariableNames", var_names);
            var_names = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
            intconv = eval(get_param(block, 'intconv_data'));
            Tintconv = array2table(intconv, "VariableNames", var_names);
            
        otherwise
            [~, ~, ~, ~, ~, Tintconv, ~, Tline] = fun_read_mask_tables(gcb);
    end
    [~, ~, ~, ~, ~, ~, Tbus1, Tline1] = fun_read_mask_tables(gcb);
    Tbus.Bs(startsWith(Tbus1.Complexity_bus, 'DC ')) = Tbus1.Bs(startsWith(Tbus1.Complexity_bus, 'DC '));
    % Tline.x(startsWith(Tline1.Complexity, 'DC ')) = Tline1.x(startsWith(Tline1.Complexity, 'DC '));
    % Tline.b(startsWith(Tline1.Complexity, 'DC ')) = Tline1.b(startsWith(Tline1.Complexity, 'DC '));
    fun_fill_tables_in_mask(gcb, Tbus, Tgen, Tline, Tintconv)
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end
