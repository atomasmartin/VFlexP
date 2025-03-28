function fun_mask_delete_row_Fb(block)
    maskObj = Simulink.Mask.get(block);
    tableControl = maskObj.getDialogControl( 'Fb_data' );
    rowIndex = tableControl.getSelectedRows();
    if ~isempty(rowIndex)
        tableControl.removeRow(rowIndex(1));
    end
end