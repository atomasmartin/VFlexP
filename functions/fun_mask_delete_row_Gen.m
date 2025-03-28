function fun_mask_delete_row_Gen(block)
    maskObj = Simulink.Mask.get(block);
    tableControl = maskObj.getDialogControl( 'gen_data' );
    rowIndex = tableControl.getSelectedRows();
    if ~isempty(rowIndex)
        tableControl.removeRow(rowIndex(1));
    end
end