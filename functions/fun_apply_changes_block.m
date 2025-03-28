function fun_apply_changes_block(gcb)
try
    selected_block_handle = get_param(gcb,'Handle');
    dlgs = DAStudio.ToolRoot.getOpenDialogs;
    for i=1:length(dlgs) %find dialog of selected block
        if class(dlgs(i).getSource) == "Simulink.SLDialogSource"
            dialog_block_handle = dlgs(i).getSource.getBlock.Handle;
            if dialog_block_handle == selected_block_handle
                dlgs(i).apply %'click' apply
            end
        end
    end
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end