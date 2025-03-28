function fun_mask_new_row_Param(block)
    maskObj = Simulink.Mask.get(block);
    tableControl = maskObj.getDialogControl('param_change_data');
    rowIndex = tableControl.getSelectedRows();
    if ~isempty(rowIndex)
        rowIndex = rowIndex(1);
        % Store the values of the selected row to copy them to the new row
        values_str = strings(1, tableControl.getNumberOfColumns);
        for i=1:length(values_str)
            values_str(i) = tableControl.getValue([rowIndex, i]);
        end
        % Add a new row
        if rowIndex ~= tableControl.getNumberOfRows
            eval(strcat("tableControl.insertRow(rowIndex+1, ", join(repmat("''", 1, tableControl.getNumberOfColumns), " , "), ")"))
        else
            eval(strcat("tableControl.addRow(", join(repmat("''", 1, tableControl.getNumberOfColumns), " , "), ")"))
        end
        for i=1:length(values_str)
            tableControl.setValue([rowIndex+1, i], values_str(i));
        end
    else
        % Add a new row
        eval(strcat("tableControl.addRow(", join(repmat("''", 1, tableControl.getNumberOfColumns), " , "), ")"))
    end
end