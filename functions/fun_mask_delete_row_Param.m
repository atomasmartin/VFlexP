function fun_mask_delete_row_Param(block)
    maskObj = Simulink.Mask.get(block);
    tableControl = maskObj.getDialogControl( 'param_change_data' );
    rowIndex = tableControl.getSelectedRows();
    if ~isempty(rowIndex)
        tableControl.removeRow(rowIndex(1));
    end
end