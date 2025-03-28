function fun_mask_delete_row_Line(block)
    maskObj = Simulink.Mask.get(block);
    tableControl = maskObj.getDialogControl( 'line_data' );
    rowIndex = tableControl.getSelectedRows();
    if ~isempty(rowIndex)
        tableControl.removeRow(rowIndex(1));
    end
end