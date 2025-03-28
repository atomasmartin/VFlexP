function fun_initialise_system_from_tables(gcb)
try
    block = gcb;
    mdlWks = get_param(bdroot, "ModelWorkspace");
    [Tbus, Tgen, Tline, TbusDC, TlineDC, Tintconv, ~, ~, TFb] = fun_read_mask_tables(block);
    bvarname = get_param(block, "bvarname");
    fun_load_complexity_grid(gcb)

    % Save gen_static_params
    gen_static_params.mBase = Tgen.mBase;
    assignin(mdlWks, [bvarname 'gen_static_params'], gen_static_params)

    % set one reference generator per each area
    fref.areas = TFb.Area;
    fref.ref_units = zeros(1, length(fref.areas));
    for i = 1:length(fref.areas)
        fref.ref_units(i) = find(ismember(Tgen.bus, Tbus.bus_i(Tbus.area == fref.areas(i))), 1);
    end

    % Calculate base frequency for every unit
    fref.Fb_units = zeros(1, height(Tgen));
    for i = 1:length(fref.Fb_units)
        fref.Fb_units(i) = find(fref.areas == Tbus.area(Tbus.bus_i == Tgen.bus(i)), 1);
    end

    fref.buses = zeros(1, height(Tbus));
    for i = 1:height(Tbus)
        fref.buses(i) = find(fref.areas == Tbus.area(i));
    end
    fref.lines = zeros(1, height(Tline));
    for i = 1:height(Tline)
        fref.lines(i) = find(fref.areas == Tbus.area(Tbus.bus_i == Tline.fbus(i)));
    end
    loads_idx = find(Tbus.Pd ~= 0 | Tbus.Qd ~= 0);
    fref.loads = zeros(1, length(loads_idx));
    for i = 1:length(loads_idx)
        fref.loads(i) = find(fref.areas == Tbus.area(loads_idx(i)));
    end
    % Units
    Gen = unique(Tgen.type);
    for i = 1:length(Gen)
        idx_Gen = find(strcmp(Tgen.type, Gen(i)));
        fref.(Gen(i)) = zeros(1, length(idx_Gen));
        for ii = 1:length(idx_Gen)
            fref.(Gen(i))(ii) = find(fref.areas == Tbus.area(Tbus.bus_i == Tgen.bus(idx_Gen(ii))), 1);
        end
    end
    % DC-AC converters
    if ~isempty(Tintconv)
        Gen = unique(Tintconv.type);
        for i = 1:length(Gen)
            idx_Gen = find(strcmp(Tintconv.type, Gen(i)));
            fref.(Gen(i)) = zeros(1, length(idx_Gen));
            for ii = 1:length(idx_Gen)
                fref.(Gen(i))(ii) = find(fref.areas == Tbus.area(Tbus.bus_i == Tintconv.tbus(idx_Gen(ii))));
            end
        end
    end
    assignin(mdlWks, [bvarname 'fref'], fref) % Assign to model workspace

    % If loads are constant PQ on PF, upload constant impedance for
    % equal initial PQ in dynamic

    Sb = evalin("base", get_param(gcb, 'Sb'));
    Load_buses = find(Tbus.Pd ~= 0 | Tbus.Qd ~= 0);
    if isempty(Load_buses)
        R_Loads = 0;
        L_Loads = 0;
        assignin(mdlWks, [bvarname 'Ploads'], [])
        assignin(mdlWks, [bvarname 'Qloads'], [])
    else
        Ploads = Tbus.Pd(Load_buses)';
        Qloads = Tbus.Qd(Load_buses)';
        Ploads(Ploads == 0 & Qloads == 0) = 1e-10;
        Qloads(Ploads == 0 & Qloads == 0) = 1e-10;
        R_Loads = Tbus.Vm(Load_buses)'.^2.*Sb.*(Ploads./(Ploads.^2+Qloads.^2));
        L_Loads = Tbus.Vm(Load_buses)'.^2.*Sb.*(Qloads./(Ploads.^2+Qloads.^2));
        assignin(mdlWks, [bvarname 'Ploads'], Ploads)
        assignin(mdlWks, [bvarname 'Qloads'], Qloads)
    end
    assignin(mdlWks, [bvarname 'Load_buses'], Load_buses);
    assignin(mdlWks, [bvarname 'R_Loads'], R_Loads);
    assignin(mdlWks, [bvarname 'L_Loads'], L_Loads);
    

    aux = zeros(height(Tgen), 1);
    for i=1:length(aux)
        if isempty(find(Tbus.bus_i == Tgen.bus(i), 1))
            errordlg(['Generator ' int2str(i) ' in non-existent bus'])
        end
        aux(i) = find(Tbus.bus_i == Tgen.bus(i));
    end
    buses_lf = aux;
    assignin(mdlWks, [bvarname 'buses_lf'], buses_lf) % Assign to model workspace

    Ps_lf = Tgen.Pg'/Sb; % Active power in pu
    Qs_lf = Tgen.Qg'/Sb; % Reactive power in pu
    Ps_lf(Ps_lf == 0) = 1e-18;
    Qs_lf(Qs_lf == 0) = 1e-18;
    iPs_lf = 1./Ps_lf; % Inverse of active power in pu
    iQs_lf = 1./Qs_lf; % Inverse of reactive power in pu
    iPs_lf(isinf(iPs_lf)) = 1e18;
    iQs_lf(isinf(iQs_lf)) = 1e18;
    assignin(mdlWks, [bvarname 'iPs_lf'], iPs_lf) % Assign to model workspace
    assignin(mdlWks, [bvarname 'iQs_lf'], iQs_lf) % Assign to model workspace
    assignin(mdlWks, [bvarname 'Ps_lf'], Ps_lf) % Assign to model workspace
    assignin(mdlWks, [bvarname 'Qs_lf'], Qs_lf) % Assign to model workspace

    if ~isempty(Tline)
        from = Tline.fbus'; % From bus
        to = Tline.tbus'; % To bus
        R_pu = Tline.r; % Resistance
        X_pu = Tline.x; % Reactance
        C_pu = Tline.b; % Capacitance
    end

    C_to_pu = Tbus.Bs/Sb; % Capacitor to ground at to bus
    PQ_gen = Tbus.type(buses_lf); % Bus type

    assignin(mdlWks, [bvarname 'DGs'], Tgen.type) % Assign to model workspace
    assignin(mdlWks, [bvarname 'bus'], Tbus.bus_i) % Assign to model workspace
    assignin(mdlWks, [bvarname 'PQ_gen'], PQ_gen) % Assign to model workspace

    v_buses_mod = Tbus.Vm'; % Voltage mod
    v_buses_angle = deg2rad(Tbus.Va)'; % Voltage angle
    assignin(mdlWks, [bvarname 'v_buses_mod'], v_buses_mod) % Assign to model workspace

    bus = Tbus.bus_i;
    Nbuses = length(Tbus.bus_i); % Count the number of buses

    DGs_type = unique(Tgen.type);

    for i=1:length(DGs_type)
        % Get the buses where units are connected
        assignin(mdlWks, [bvarname char(DGs_type(i)) 's_buses'], buses_lf(strcmp(Tgen.type, DGs_type(i))))
    end


    assignin(mdlWks,[bvarname 'Nbuses'], Nbuses) % Assign to the model workspace

    if ~isempty(Tline)
        aux = from;
        for i=1:length(from)
            from(i) = find(bus == aux(i), 1);
        end
        assignin(mdlWks, [bvarname 'From_buses'], from)
        aux = to;
        for i=1:length(to)
            to(i) = find(bus == aux(i), 1);
        end
        assignin(mdlWks, [bvarname 'To_buses'], to)

        Lines_toM = zeros(Nbuses, length(from));
        Lines_fromM = Lines_toM;
        for i=1:length(from)
            Lines_toM(to(i), i) = 1;
            Lines_fromM(from(i), i) = 1;
        end
        assignin(mdlWks, [bvarname 'Lines_fromM'], sparse(Lines_fromM))
        assignin(mdlWks, [bvarname 'Lines_toM'], sparse(Lines_toM))

        assignin(mdlWks, [bvarname 'R_Lines'], R_pu')
        assignin(mdlWks, [bvarname 'L_Lines'], X_pu')
        

        tap_lines = Tline.ratio;
        tap_lines(tap_lines == 0) = 1;
        shift_lines = deg2rad(Tline.angle);
        assignin(mdlWks, [bvarname 'tap_lines'], tap_lines')
        assignin(mdlWks, [bvarname 'shift_lines'], shift_lines')
    end

    assignin(mdlWks, [bvarname 'Rcf'], 1./Tbus.Gs'*Sb)
    % Cbuses:
    Cbuses = zeros(1, Nbuses);
    if ~isempty(Tline)
        for i=1:Nbuses
            Cbuses(i) = sum(C_pu(to == i | from == i))/2 + C_to_pu(i);
        end
    else
        Cbuses = C_to_pu;
    end
    assignin(mdlWks, [bvarname 'Cbuses'], Cbuses)
    assignin(mdlWks, [bvarname 'C_Lines'], C_pu');

    lf_op.op_report.Outputs.y = Ps_lf; % P gen
    lf_op.op_report.Outputs(2).y = Qs_lf; % Q gen

    if ~isempty(Tline)
        vfrom = v_buses_mod(from).*exp(v_buses_angle(from)*1i);
        vfrom = vfrom./tap_lines'.*exp(-shift_lines'.*1i); % including transformer
        v_lines_d = real(vfrom-v_buses_mod(to).*exp(v_buses_angle(to)*1i));
        v_lines_q = imag(vfrom-v_buses_mod(to).*exp(v_buses_angle(to)*1i));
        i_lines_d = (1./(R_pu'.^2+X_pu'.^2)).*(R_pu'.*v_lines_d + X_pu'.*v_lines_q);
        i_lines_q = (1./(R_pu'.^2+X_pu'.^2)).*(R_pu'.*v_lines_q - X_pu'.*v_lines_d);
        lf_op.op_report.Outputs(3).y = i_lines_d;
        lf_op.op_report.Outputs(4).y = i_lines_q;
    else
        lf_op.op_report.Outputs(3).y = [];
        lf_op.op_report.Outputs(4).y = [];
    end
    if ~isempty(Load_buses)
        v_loads_d = real(v_buses_mod(Load_buses).*exp(v_buses_angle(Load_buses)*1i));
        v_loads_q = imag(v_buses_mod(Load_buses).*exp(v_buses_angle(Load_buses)*1i));
        i_loads_d = (1./(R_Loads.^2+L_Loads.^2)).*(R_Loads.*v_loads_d + L_Loads.*v_loads_q);
        i_loads_q = (1./(R_Loads.^2+L_Loads.^2)).*(R_Loads.*v_loads_q - L_Loads.*v_loads_d);
        lf_op.op_report.Outputs(5).y = i_loads_d;
        lf_op.op_report.Outputs(6).y = i_loads_q;
        lf_op.op_report.Outputs(7).y = real(v_buses_mod.*exp(v_buses_angle*1i));
        lf_op.op_report.Outputs(8).y = imag(v_buses_mod.*exp(v_buses_angle*1i));
    else
        lf_op.op_report.Outputs(5).y = real(v_buses_mod.*exp(v_buses_angle*1i));
        lf_op.op_report.Outputs(6).y = imag(v_buses_mod.*exp(v_buses_angle*1i));
    end

    assignin(mdlWks, [bvarname 'lf_op'], lf_op)

    set_param([block '/Dynamic/Grid/Dynamic/Lines/dynamic/i_d'], 'InitialCondition', ...
        strcat(bvarname, "lf_op.op_report.Outputs(3).y(line_complexity == 1)")) % Initialise lines
    set_param([block '/Dynamic/Grid/Dynamic/Lines/dynamic/i_q'], 'InitialCondition', ...
        strcat(bvarname, "lf_op.op_report.Outputs(4).y(line_complexity == 1)")) % Initialise lines
    if ~isempty(Load_buses)
        set_param([block '/Dynamic/Grid/Dynamic/Loads/series RL dynamic/i_d'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(5).y(load_complexity == 1)")) % Initialise loads
        set_param([block '/Dynamic/Grid/Dynamic/Loads/series RL dynamic/i_q'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(6).y(load_complexity == 1)")) % Initialise loads

        set_param([block '/Dynamic/Grid/Dynamic/Loads/parallel RL dynamic/i_d'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(8).y(Load_buses(load_complexity == 6)).*Qloads" + ...
            "(load_complexity == 6)./(v_buses_mod(Load_buses(load_complexity == 6)).^2)/Sb")) % Initialise loads
        set_param([block '/Dynamic/Grid/Dynamic/Loads/parallel RL dynamic/i_q'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(7).y(Load_buses(load_complexity == 6)).*-Qloads" + ...
            "(load_complexity == 6)./(v_buses_mod(Load_buses(load_complexity == 6)).^2)/Sb")) % Initialise loads


        set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic PQ/i_d/I'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(5).y(load_complexity == 5)")) % Initialise loads
        set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic PQ/i_q/I'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(6).y(load_complexity == 5)")) % Initialise loads
        set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vd'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(7).y(bus_complexity == 1)")) % Initialise buses
        set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vq'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(8).y(bus_complexity == 1)")) % Initialise buses
    else
        set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vd'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(5).y(bus_complexity == 1)")) % Initialise buses
        set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vq'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(6).y(bus_complexity == 1)")) % Initialise buses
    end

    if ~isempty(TbusDC) % DC grid
        DC.include = true;
        DC.Nbuses = height(TbusDC); % Assign to the model workspace
        DC.from = TlineDC.fbus;
        DC.to = TlineDC.tbus;
        DC.bus = TbusDC.bus_i;
        DC.v = TbusDC.Vm';
        [G_dc, Gnoloads_dc] = fun_admitt(block, TbusDC, TlineDC);
        
        assignin(mdlWks, [bvarname 'conductance_matrix'], real(G_dc)) % Assign to model workspace
        if rank(real(G_dc)) < length(real(G_dc))
            assignin(mdlWks, [bvarname 'iconductance_matrix'], pinv(real(G_dc))) % Assign to model workspace
        else
            assignin(mdlWks, [bvarname 'iconductance_matrix'], inv(real(G_dc))) % Assign to model workspace
        end
        assignin(mdlWks, [bvarname 'conductance_matrix_noloads'], real(Gnoloads_dc)); % Assign to model workspace
        assignin(mdlWks, [bvarname 'iconductance_matrix_noloads'], inv(real(Gnoloads_dc))) % Assign to model workspace

        Load_busesDC = find(TbusDC.Pd ~= 0);
        if isempty(Load_busesDC)
            R_Loads_DC = 0;
            assignin(mdlWks, [bvarname 'PloadsDC'], [])
        else
            PloadsDC = TbusDC.Pd(Load_busesDC)';
            PloadsDC(PloadsDC == 0) = 1e-10;
            R_Loads_DC = TbusDC.Vm(Load_busesDC)'.^2.*Sb.*(PloadsDC./(PloadsDC.^2));
            assignin(mdlWks, [bvarname 'PloadsDC'], PloadsDC)
        end
        assignin(mdlWks, [bvarname 'Load_busesDC'], Load_busesDC);
        assignin(mdlWks, [bvarname 'R_Loads_DC'], R_Loads_DC);

        aux = DC.from;
        for i=1:length(DC.from)
            DC.from(i) = find(DC.bus == aux(i), 1);
        end
        aux = DC.to;
        for i=1:length(DC.to)
            DC.to(i) = find(DC.bus == aux(i), 1);
        end
        Lines_toM = zeros(DC.Nbuses, length(DC.from));
        Lines_fromM = Lines_toM;
        for i=1:length(DC.from)
            Lines_toM(DC.to(i), i) = 1;
            Lines_fromM(DC.from(i), i) = 1;
        end
        DC.Lines_toM = sparse(Lines_toM);
        DC.Lines_fromM = sparse(Lines_fromM);
        if ~isempty(TlineDC)
            DC.R_Lines = TlineDC.r';
            DC.L_Lines = TlineDC.x';
            DC.CbusDC = TbusDC.Bs;
            NbusesDC = height(TbusDC);
            for i=1:NbusesDC
                DC.CbusDC(i) = sum(TlineDC.b(DC.to == i | DC.from == i))/2 + TbusDC.Bs(i);
            end
        else
            DC.R_Lines = [];
        end
        [~, indices_from] = ismember(TlineDC.fbus, TbusDC.bus_i);
        [~, indices_to] = ismember(TlineDC.tbus, TbusDC.bus_i);
        DC.i_branch = (TbusDC.Vm(indices_from) - TbusDC.Vm(indices_to))./DC.R_Lines;
        DC.Rbus = Sb./TbusDC.Gs;
    else
        DC.include = false;
        intconv_PF.include = false;
        assignin(mdlWks, [bvarname 'intconv_PF'], intconv_PF)
    end
    assignin(mdlWks, [bvarname 'DC'], DC)
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end
