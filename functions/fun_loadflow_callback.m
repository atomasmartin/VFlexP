function fun_loadflow_callback(gcb)
if ~bdIsLibrary(bdroot)
    try

        wb = waitbar(0, 'PF', 'Name', ...
            'Initialising the simulink file'); % Waitbar
        block = gcb; % Get the current block name
        bvarname = get_param(block, 'bvarname');
        [Tbus, ~, Tline] = fun_read_mask_tables(block);
        [y, ynoloads] = fun_admitt(block, Tbus, Tline);
        model = bdroot;
        mdlWks = get_param(model,'ModelWorkspace'); % Get the model workspace
        assignin(mdlWks, [bvarname 'admitt_matrix'], y) % Assign to model workspace
        assignin(mdlWks, [bvarname 'iadmitt_matrix'], inv(y)) % Assign to model workspace
        assignin(mdlWks, [bvarname 'admitt_matrix_noloads'], ynoloads) % Assign to model workspace
        assignin(mdlWks, [bvarname 'iadmitt_matrix_noloads'], inv(ynoloads)) % Assign to model workspace
        
        set_param(block, 'LabelModeActiveChoice', 'Loadflow') % Set active choice to Loadflow
        close_system(block)
        open_system(block, 'mask')
        [Tbus, Tgen, Tline] = fun_read_mask_tables(gcb);
        fun_set_complexity_blocks_grid(gcb, Tline, Tbus)
        if height(Tbus) == 0 || height(Tgen) == 0 || height(Tline) == 0
            return
        end
        waitbar(0.1,wb)

        Sb = evalin("base", get_param(gcb, 'Sb'));

        Load_buses = find(Tbus.Pd ~= 0 | Tbus.Qd ~= 0);
        noLoads = isempty(Load_buses);
        if noLoads
            R_Loads = 0;
            L_Loads = 0;
            assignin(mdlWks, [bvarname 'Ploads'], [])
            assignin(mdlWks, [bvarname 'Qloads'], [])
        else
            Ploads = Tbus.Pd(Load_buses)';
            Qloads = Tbus.Qd(Load_buses)';
            Ploads(Ploads == 0) = 1;
            Qloads(Qloads == 0) = 1;
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
            aux(i) = find(Tbus.bus_i == Tgen.bus(i));
            if isempty(aux(i))
                errordlg(['Generator ' int2str(i) ' in non-existent bus'])
            end
        end
        buses_lf = aux;
        assignin(mdlWks, [bvarname 'buses_lf'], buses_lf) % Assign to model workspace

        Ps_lf = Tgen.Pg'/Sb; % Active power in pu
        Qs_lf = Tgen.Qg'/Sb; % Reactive power in pu
        Ps_lf(Ps_lf == 0) = 1e-6;
        Qs_lf(Qs_lf == 0) = 1e-6;
        iPs_lf = 1./Ps_lf; % Inverse of active power in pu
        iQs_lf = 1./Qs_lf; % Inverse of reactive power in pu
        iPs_lf(isinf(iPs_lf)) = 1e18;
        iQs_lf(isinf(iQs_lf)) = 1e18;
        assignin(mdlWks, [bvarname 'iPs_lf'], iPs_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'iQs_lf'], iQs_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'Ps_lf'], Ps_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'Qs_lf'], Qs_lf) % Assign to model workspace

        from = Tline.fbus'; % From bus
        to = Tline.tbus'; % To bus
        R_pu = Tline.r; % Resistance
        X_pu = Tline.x; % Reactance
        C_pu = Tline.b; % Capacitance

        C_to_pu = Tbus.Bs/Sb; % Capacitor to ground at to bus
        PQ_gen = Tbus.type(buses_lf)'; % Bus type

        assignin(mdlWks, [bvarname 'DGs'], Tgen.type) % Assign to model workspace
        assignin(mdlWks, [bvarname 'bus'], Tbus.bus_i) % Assign to model workspace
        assignin(mdlWks, [bvarname 'PQ_gen'], PQ_gen) % Assign to model workspace

        v_buses_mod = Tbus.Vm'; % Voltage mod
        assignin(mdlWks, [bvarname 'v_buses_mod'], v_buses_mod) % Assign to model workspace

        bus = Tbus.bus_i;
        Nbuses = length(Tbus.bus_i); % Count the number of buses

        DGs_type = unique(Tgen.type);

        for i=1:length(DGs_type)
            % Get the buses where units are connected
            assignin(mdlWks, [bvarname char(DGs_type(i)) 's_buses'], buses_lf(strcmp(Tgen.type, DGs_type(i))))
        end

        assignin(mdlWks,[bvarname 'Nbuses'], Nbuses) % Assign to the model workspace
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
        assignin(mdlWks, [bvarname 'R_Lines'], R_pu')
        assignin(mdlWks, [bvarname 'L_Lines'], X_pu')

        tap_lines = Tline.ratio;
        tap_lines(tap_lines == 0) = 1;
        shift_lines = deg2rad(Tline.angle);
        assignin(mdlWks, [bvarname 'tap_lines'], tap_lines')
        assignin(mdlWks, [bvarname 'shift_lines'], shift_lines')

        assignin(mdlWks, [bvarname 'Rcf'], 1./Tbus.Gs'*Sb)
        % Cbuses:
        Cbuses = zeros(1, Nbuses);
        for i=1:Nbuses
            Cbuses(i) = sum(C_pu(to == i | from == i))/2 + C_to_pu(i);
        end
        assignin(mdlWks, [bvarname 'Cbuses'], Cbuses)


        waitbar(0.3,wb)

        opspec = operspec(model); % Create operating point specification
        opspec = addoutputspec(opspec, [block '/Loadflow/Generators/PC'], 3); % P
        waitbar(0.4,wb)
        opspec = addoutputspec(opspec, [block '/Loadflow/Generators/PC'], 4); % Q
        opspec = addoutputspec(opspec, [block '/Loadflow/Generators/vd'], 1); % Buses vd
        opspec = addoutputspec(opspec, [block '/Loadflow/Generators/vq'], 1); % Buses vq

        waitbar(0.7,wb)

        optionop = findopOptions('OptimizerType', 'graddescent-elim', 'OptimizationOptions', optimset('PlotFcns', ...
            'optimplotfirstorderopt', "MaxFunEvals", ...
            2000*Nbuses, "MaxIter", 50));
        optionop.OptimizationOptions.Algorithm = 'interior-point';
        [lf_op.op, lf_op.op_report] = findop(model, opspec, optionop); % Find operating point
        close % close PlotFcns
        if ~strcmp(lf_op.op_report.TerminationString, "Operating point specifications were successfully met.")
            close(wb)
            assignin(mdlWks, [bvarname 'lf_op'], lf_op) % Assign to model workspace
            errordlg("Generate loadflow report for details", "Loadflow did not converge")
            return
        end
        % If loads are constant PQ on PF, upload constant impedance for
        % equal initial PQ in dynamic
        if noLoads == 1
            R_Loads = 0;
            L_Loads = 0;
        else
            vmod = abs(lf_op.op_report.Outputs(3).y(Load_buses) + ...
                lf_op.op_report.Outputs(4).y(Load_buses)*1i)';
            R_Loads = vmod.^2.*Sb.*(Ploads./(Ploads.^2+Qloads.^2));
            L_Loads = vmod.^2.*Sb.*(Qloads./(Ploads.^2+Qloads.^2));
        end
        assignin(mdlWks, [bvarname 'R_Loads'], R_Loads); % Convert to pu
        assignin(mdlWks, [bvarname 'L_Loads'], L_Loads); % Convert to pu

        Ps_lf = lf_op.op_report.Outputs(1).y'; % Active power in pu
        Qs_lf = lf_op.op_report.Outputs(2).y'; % Reactive power in pu
        iPs_lf = 1./Ps_lf; % Inverse of active power in pu
        iQs_lf = 1./Qs_lf; % Inverse of reactive power in pu
        iPs_lf(isinf(iPs_lf)) = 1e18;
        iQs_lf(isinf(iQs_lf)) = 1e18;
        assignin(mdlWks, [bvarname 'iPs_lf'], iPs_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'iQs_lf'], iQs_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'Ps_lf'], Ps_lf) % Assign to model workspace
        assignin(mdlWks, [bvarname 'Qs_lf'], Qs_lf) % Assign to model workspace

        set_param([block '/Dynamic/Grid/Dynamic/Lines/dynamic/i_d'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(3).y(line_complexity == 1)", "'")) % Initialise lines
        set_param([block '/Dynamic/Grid/Dynamic/Lines/dynamic/i_q'], 'InitialCondition', ...
            strcat(bvarname, "lf_op.op_report.Outputs(4).y(line_complexity == 1)", "'")) % Initialise lines
        if noLoads == 0
            set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic/i_d'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(5).y(load_complexity == 1)", "'")) % Initialise loads
            set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic/i_q'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(6).y(load_complexity == 1)", "'")) % Initialise loads
            set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic PQ/i_d/I'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(5).y(load_complexity == 5)", "'")) % Initialise loads
            set_param([block '/Dynamic/Grid/Dynamic/Loads/dynamic PQ/i_q/I'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(6).y(load_complexity == 5)", "'")) % Initialise loads
            set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vd'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(7).y(bus_complexity == 1)", "'")) % Initialise buses
            set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vq'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(8).y(bus_complexity == 1)", "'")) % Initialise buses
        else
            set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vd'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(5).y(bus_complexity == 1)", "'")) % Initialise buses
            set_param([block '/Dynamic/Grid/Dynamic/Buses/dynamic/vq'], 'InitialCondition', ...
                strcat(bvarname, "lf_op.op_report.Outputs(6).y(bus_complexity == 1)", "'")) % Initialise buses
        end


        assignin(mdlWks, [bvarname 'lf_op'], lf_op) % Assign to model workspace

        set_param(block, 'LabelModeActiveChoice', 'Dynamic') % Set active choice to grid
        close(wb) % Close waitbar

        close_system(block)
        % Upload table with LF solution
        v_buses = lf_op.op_report.Outputs(3).y + 1i*lf_op.op_report.Outputs(4).y;

        Tbus.Vm = abs(v_buses);
        Tbus.Va = rad2deg(angle(v_buses));
        slack_idx = find(Tbus.type == 3, 1);
        if isempty(slack_idx)
            slack_idx = 1;
        end
        Tbus.Va = wrapTo180(Tbus.Va - Tbus.Va(slack_idx));
        Tgen.Pg = Sb.*lf_op.op_report.Outputs(1).y;
        Tgen.Qg = Sb.*lf_op.op_report.Outputs(2).y;
        fun_fill_tables_in_mask(block, Tbus, Tgen, Tline, [])
        fun_initialise_system_from_tables(block)
        open_system(block, 'mask')

    catch ME
        str_name = string({ME.stack.name}');
        str_line = string({ME.stack.line}');
        errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
        errordlg(char(errortext));
    end
end
end