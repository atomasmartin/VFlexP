function fun_plot_change_parameter_callback(gcb)
try
model = bdroot;
block = gcb; % Get the current block name

mdlWks = get_param(model,'ModelWorkspace'); % Get the model workspace

maskObj = Simulink.Mask.get(gcb);
tableControl = maskObj.getDialogControl( 'param_change_data' );
rowIndex = tableControl.getSelectedRows();

param_change_data_table = eval(get_param(block, 'param_change_data')); % Get the param change table


if ~isempty(rowIndex)
    
    Real_Z = [];
    Imag_Z = [];
    LinsysParamChangeStruct = evalin(mdlWks,'LinsysParamChangeStruct');
    FromValueToChange = eval(cell2mat(param_change_data_table(rowIndex,2)));
    ToValueToChange = eval(cell2mat(param_change_data_table(rowIndex,3)));
    StepsToChange = eval(cell2mat(param_change_data_table(rowIndex,4)));
    ValuesVector = linspace(FromValueToChange,rowIndex,StepsToChange);
    ParameterSymbol = char(param_change_data_table(rowIndex,1));
    length_linsys = length(fieldnames(LinsysParamChangeStruct.("Param"+rowIndex)))

    for i = 1:numel(fieldnames(LinsysParamChangeStruct.("Param"+rowIndex)))
        StepVariable_names = fieldnames(LinsysParamChangeStruct.("Param"+rowIndex))
        StepVar = StepVariable_names(i)
        Real_Z = [Real_Z real(eig(LinsysParamChangeStruct.("Param"+rowIndex).(string(StepVar)).A))];
        Imag_Z = [Imag_Z imag(eig(LinsysParamChangeStruct.("Param"+rowIndex).(string(StepVar)).A))];
    end
    color_vector = hsv(length_linsys);

    figure(1);clf;
    hold on

    for i = 1:numel(fieldnames(LinsysParamChangeStruct.("Param"+rowIndex)))
        if i == 1
            p1 = plot(Real_Z(:,i),Imag_Z(:,i),'o','color',color_vector(i,:));
        else
            plot(Real_Z(:,i),Imag_Z(:,i),'x','color',color_vector(i,:))
        end
    end
    legend({'INITAL POINT'})
    cb = colorbar;
    col = colormap(color_vector);
    cb.FontSize = 12;        
    cb.Ticks = linspace(0,1,size(Real_Z,2)/10);
    cb.TickLabels =  linspace(FromValueToChange,ToValueToChange,size(Real_Z,2)/10);
    cb.Label.String =ParameterSymbol;       
    cb.Label.FontSize = 12;
    cb.Label.Rotation = 90;


end
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end