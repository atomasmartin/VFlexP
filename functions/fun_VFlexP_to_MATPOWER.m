function mpc = fun_VFlexP_to_MATPOWER(gcb)
    block = gcb;
    [~, ~, ~, ~, ~, Tintconv, Tbus, Tline, ~, ~, Tgen] = fun_read_mask_tables(block);
    Tbus.Bs(startsWith(Tbus.Complexity_bus, 'DC ')) = 0;
    Tline.x(startsWith(Tline.Complexity, 'DC ')) = 0;
    Tline.b(startsWith(Tline.Complexity, 'DC ')) = 0;
    if isempty(Tline)
        mpc.branch = [];
    else
        Tline.Complexity = [];
        Tline.Complexity_LF = [];
        if ~isempty(Tintconv)
            Tintconv.type = [];
            fld_names = Tintconv.Properties.VariableNames(width(Tline)+1:end);
            for i = 1:length(fld_names)
                Tline.(fld_names{i})(:) = 0;
            end
            Tline.TAP_MIN(:) = 1;
            Tline.TAP_MAX(:) = 1;
            Tline.K2(:) = 1;
            Tline.SH_MIN(:) = -360;
            Tline.SH_MAX(:) = 360;
            Tlinecomplete = [Tline; Tintconv];
        else
            Tlinecomplete = Tline;
        end
        mpc.branch = table2array(Tlinecomplete);
    end
    mpc.version = '2';
    mpc.baseMVA = evalin("base", get_param(block, 'Sb'));
    mpc.gencost = [];
    mpc.bus = table2array(Tbus(:,1:end-4));
    mpc.gen = table2array(Tgen(:,[1 3:end]));
    mpc.gencost = repmat([2	0 0 3 0 20 0], height(Tgen));
end