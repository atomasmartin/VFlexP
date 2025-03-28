function fun_load_system_callback(gcb, file)

if ~bdIsLibrary(bdroot)
    if nargin == 1
        [file, path] = uigetfile('*.xlsx','Select table file');
        full_filename = [path file];
    else
        full_filename = file;
    end

    if file ~= 0
        try
            [path,filename,ext] = fileparts(full_filename);
            % Clear all tables
            table_name = ["gen_data", "line_data", "bus_data", "intconv_data", "param_change_data"];
            for itable = 1:length(table_name)
                maskObj = Simulink.Mask.get(gcb);
                tableControl = maskObj.getDialogControl(table_name(itable));
                for i=1:tableControl.getNumberOfRows()
                    tableControl.removeRow(1);
                end
            end
            fun_apply_changes_block(gcb)
            if strcmp(ext, ".m") % MATPOWER file
                % change to that directory
                old_path = pwd;
                cd(path)
                mpc = feval(filename);
                cd(old_path)
                % get MATPOWER file
                fun_MATPOWER_to_VFlexP(mpc, gcb, 1)

            else % Excel file
                % change to that directory
                old_path = pwd;
                cd(path)
                sh_names = sheetnames(file);
                Tbus = readtable(file, "Sheet", "Bus data");
                if ismember("Line data", sh_names)
                    Tline = readtable(file, "Sheet", "Line data");
                    if isnumeric(Tline.Complexity)
                        Tline.Complexity = strings(height(Tline), 1);
                        Tline.Complexity_LF = Tline.Complexity;
                    end
                else
                    Tline = [];
                end
                Tgen = readtable(file, "Sheet", "Gen data");
                if ismember("IntConv data", sh_names)
                    Tintconv = readtable(file, "Sheet", "IntConv data");
                else
                    Tintconv = [];
                end
                if ismember("Parameter variation", sh_names)
                    Tparam = readtable(file, "Sheet", "Parameter variation");
                else
                    Tparam = [];
                end
                if ismember("System parameters", sh_names)
                    Tsysparam = readtable(file, "Sheet", "System parameters");
                    try
                        set_param(gcb, 'Sb', num2str(evalin('base', num2str((Tsysparam.Sb{1})))))
                    catch
                        set_param(gcb, 'Sb', num2str(evalin('base', num2str((Tsysparam.Sb(1))))))
                    end
                    Tsysparam = Tsysparam(:,["Area", "Fb"]);
                    % Bus table
                    str = '{';
                    Rows = height(Tsysparam);
                    Columns = width(Tsysparam);
                    for i = 1:Rows
                        for ii = 1:Columns
                            str = [str char(" '") char(string(Tsysparam{i,ii})) char("' ")];
                        end
                        if i < Rows
                            str = [str,';'];
                        end
                    end
                    str = [str,'}'];
                    set_param(gcb, 'Fb_data', str)
                end
                if ismember("Secondary Control", sh_names)
                    Tsec = readtable(file, "Sheet", "Secondary Control");
                    set_param(gcb, 'SC_structure', Tsec.SC_structure)
                    set_param(gcb, 'define_SC_in_MATLAB', evalin('base', num2str(Tsec.Define_SC_in_MATLAB{1})))
                end
                if ismember("Tertiary Control", sh_names)
                    Tter = readtable(file, "Sheet", "Tertiary Control");
                    set_param(gcb, 'TC_structure', Tsec.TC_structure)
                end
                
                if isnumeric(Tbus.Complexity_bus)
                    Tbus.Complexity_bus = strings(height(Tbus), 1);
                    Tbus.Complexity_bus_LF = Tbus.Complexity_bus;
                    Tbus.Complexity_load = Tbus.Complexity_bus;
                    Tbus.Complexity_load_LF = Tbus.Complexity_bus;
                end

                cd(old_path)
                fun_fill_tables_in_mask(gcb, Tbus, Tgen, Tline, Tintconv, Tparam)
            end

        catch ME
            str_name = string({ME.stack.name}');
            str_line = string({ME.stack.line}');
            errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
            errordlg(char(errortext));
        end
    end
end
end