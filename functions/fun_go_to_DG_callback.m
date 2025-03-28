function fun_go_to_DG_callback(gcb)
    maskObj = Simulink.Mask.get(gcb);
    tableControl = maskObj.getDialogControl( 'gen_data' );
    rowIndex = tableControl.getSelectedRows();
    grid_data_table = eval(get_param(gcb, 'gen_data'));
    DGi = grid_data_table(rowIndex, strcmp...
        ({tableControl.Columns.Name}, 'type'));
    if strlength(DGi) > 0
        if startsWith(DGi, 'DC')
            open_system([gcb '/Dynamic/Generators_DC/' DGi{1}],"tab")
        else
            open_system([gcb '/Dynamic/Generators/' DGi{1}],"tab")
            
        end
    else
        errordlg('Please select a row with a generator type', 'Select valid row')
    end
end
    