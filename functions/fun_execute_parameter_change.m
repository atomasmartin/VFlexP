function fun_execute_parameter_change(gcb)
try
    if ~bdIsLibrary(bdroot)

        tic
        model = bdroot;
        block = gcb; % Get the current block name
        bvarname = get_param(block, 'bvarname');
        mdlWks = get_param(model,'ModelWorkspace'); % Get the model workspace

        param_change_data_table = eval(get_param(block, 'param_change_data')); % Get the param change table
        save_in_workspace_flag = strcmp(get_param(block, 'save_in_workspace'),'on');
        parameter_change_filename = ['Parameter_Change_V_' get_param(block, 'parameter_change_filename')];
        NumberOfParamsToChange = size(param_change_data_table,1);
        maskObj = Simulink.Mask.get(gcb);
        LinsysParamChangeStruct = struct();
        operspecstruct = struct();
        steps_counter = 0;
        fun_close_all_dialogs;

        % Obtener el tamaño de la pantalla
        screenSize = get(0, 'ScreenSize');  % Devuelve [left, bottom, width, height]

        for ParamItem = 1:NumberOfParamsToChange %for each parameter that needs to be changed:

            ParameterSymbolALL = string(param_change_data_table(:,1));

            for i = 1:length(string(ParameterSymbolALL))
                assignin('base',ParameterSymbolALL(i),eval(cell2mat(param_change_data_table(i,2))));
            end

            FromValueToChange = eval(cell2mat(param_change_data_table(ParamItem,3)));
            ToValueToChange = eval(cell2mat(param_change_data_table(ParamItem,4)));
            StepsToChange = eval(cell2mat(param_change_data_table(ParamItem,5)));
            ValuesVector = linspace(FromValueToChange,ToValueToChange,StepsToChange);
            ParameterSymbol = char(param_change_data_table(ParamItem,1));

            io = getlinio(model);

            for StepNumber = 1:StepsToChange
                K = ValuesVector(StepNumber);
                assignin('base',ParameterSymbol,K);


                fun_apply_changes_block(gcb)
                close_system(gcb)
                open_system(gcb, 'mask')

                fun_MATPOWER_loadflow(gcb, 0)

                if isempty(allchild(0))
                    all_ok = 1;
                else
                    if all(strcmp({allchild(0).Name},'EXECUTION TIME'))
                        all_ok = 1;
                    else
                        dialogs = allchild(0);
                        close(dialogs(~strcmp({dialogs.Name},'EXECUTION TIME')));
                        all_ok = 0;

                    end
                end
                if all_ok
                    opsmodel = operspec(model);

                    operspecstruct.("Param"+ParamItem).("Step"+StepNumber) =opsmodel;

                    linsys = linearize(model,io,opsmodel);

                    if isempty(linsys)
                        errordlg("Empty linsys")
                        return
                    end
                    if height(linsys)>1
                        linsys = linsys(1,1);
                    end
                    % for i = 1:length(string(ParameterSymbolALL))
                    %     assignin('base', ParameterSymbolALL(i),eval(cell2mat(param_change_data_table(i,2))));
                    % end

                    LinsysParamChangeStruct.("Param"+ParamItem).("Step"+StepNumber) = linsys;
                    steps_counter = steps_counter+1;

                    if ParamItem == 1 && StepNumber == 1
                        Number_of_total_steps = 0;
                        first_elapse_time = toc;
                        for k = 1:NumberOfParamsToChange
                            Number_of_total_steps = Number_of_total_steps+eval(cell2mat(param_change_data_table(k,end)));
                        end
                        Expected_Time = first_elapse_time*Number_of_total_steps*1.2; %1.2 is 20% more time due to computations
                        wb_time = waitbar(1-(Expected_Time-steps_counter*first_elapse_time)/Expected_Time,['REMAINING EXECUTION TIME:  ', char(minutes((Expected_Time-steps_counter*first_elapse_time)/60))], 'Name', ...
                            'EXECUTION TIME'); % Waitbar
                    end
                    try % this try is just in case someone closed the wait bar...
                        if ~isempty(wb_time)
                            close(wb_time);
                        end

                        wb_time = waitbar(1-(Expected_Time-steps_counter*first_elapse_time)/Expected_Time,['REMAINING EXECUTION TIME:  ', char(minutes((Expected_Time-steps_counter*first_elapse_time)/60))], 'Name', ...
                            'EXECUTION TIME'); % Waitbar
                        % Calcular la posición centrada horizontalmente y 1/3 desde la parte superior
                        msgBoxSize = get(wb_time, 'Position');  % Devuelve [left, bottom, width, height] del msgbox
                        newLeft = (screenSize(3) - msgBoxSize(3)) / 2;  % Centrar horizontalmente
                        newBottom = screenSize(4) * 2/3 - msgBoxSize(4) / 2;  % 1/3 desde la parte superior

                        % Establecer la nueva posición
                        set(wb_time, 'Position', [newLeft, newBottom, msgBoxSize(3), msgBoxSize(4)]);
                    catch
                    end

                else
                    warningDialogHandles = allchild(0);
                    close(warningDialogHandles);
                end
            end
        end
        assignin(mdlWks,'LinsysParamChangeStruct',LinsysParamChangeStruct);
        assignin('base','LinsysParamChangeStruct',LinsysParamChangeStruct);
        assignin(mdlWks,'operspecstruct',operspecstruct);
        assignin('base','operspecstruct',operspecstruct);
        elapse_param_change = toc;
        try % this try is just in case someone closed the wait bar...
            if ~isempty(wb_time)
                close(wb_time);
            end
            msg = ['Parameter variation finalized in ',char((minutes(elapse_param_change/60)))];
            h = msgbox(msg, 'Success ', 'help');

        catch
        end

        if save_in_workspace_flag

            command_1 = ['save ' parameter_change_filename];
            % command_1 = ['save '  [filename '_PolosSalientesDFIG_PQ_SC_TrafoDFIG_variation_' int2str(100*(zeta_pll))]];
            evalin("base",command_1)
        end
    end
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end