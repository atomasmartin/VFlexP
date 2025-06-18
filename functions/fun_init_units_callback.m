function fun_init_units_callback(gcb)
    if ~bdIsLibrary(bdroot)
        try

            block = gcb; % Get the current block name
            [Tbus, Tgen, Tline, ~, ~, Tintconv, ~, ~, TFb] = fun_read_mask_tables(block);
            
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
            

        wb = waitbar(0, 'Initialising', 'Name', ...
            'Initialising the units'); % Waitbar
        block = gcb; % Get the current block name
        bvarname = get_param(block, 'bvarname');
        waitbar(0.1,wb);

        [Tbus, Tgen, ~, ~, ~, ~, ~, ~, TFb, TgenDC] = fun_read_mask_tables(block);
        DGs = Tgen.type;
        Generators = unique(DGs);
        omega_base = 2*pi*TFb.("Fb")';
        mdlWks = get_param(bdroot,'ModelWorkspace'); % Get the model workspace
        assignin(mdlWks, [bvarname 'init_OK_bool'], true)
        assignin(mdlWks, [bvarname 'fref'], fref) % Assign to model workspace
        bl = find_system([block '/Dynamic/Generators'], ...
            'IncludeCommented', 'on', 'SearchDepth', '2', ...
            'LookUnderMasks', 'on', 'Name', 'trim_output');
        for i = 1:length(bl)
            bl{i} = get_param(bl{i}, 'Parent');
        end
        for i=1:length(bl) 
            if sum(strcmp(DGs, get_param(bl{i}, 'Name'))) == 0
                set_param(bl{i}, "commented", "on")
            else
                set_param(bl{i}, "commented", "off")
            end
        end   
        waitbar(0.2,wb);
        % Initialize values to guarantee compilation
        for i = 1:length(bl) 
           DG_namei = get_param(bl{i}, 'Name');
           assignin(mdlWks, [bvarname DG_namei 's'], 1)
           assignin(mdlWks, [bvarname DG_namei 's_buses'],1)
        end
        waitbar(0.4,wb);
        for i=1:length(Generators)
            % Get the buses where units are connected
            DGs_typei_location = find(strcmp(DGs, Generators{i}));
            aux = zeros(1, length(DGs_typei_location));
            for ii = 1:length(aux)
                aux(ii) = find(Tbus.bus_i == Tgen.bus(DGs_typei_location(ii)));
            end
            assignin(mdlWks, [bvarname Generators{i} 's_buses'], aux')

            aux_M = zeros(height(Tbus), length(DGs_typei_location));
            aux_idx = sub2ind(size(aux_M), aux, 1:length(DGs_typei_location));
            aux_M(aux_idx) = 1;
            assignin(mdlWks, [bvarname Generators{i} 's_toM'], aux_M)
         
        end

        waitbar(0.5,wb);
    
        for i=1:length(Generators)

            assignin(mdlWks, [bvarname Generators{i} '_idx'], find(strcmp(DGs, Generators{i}))); % Indices of unit
    
            % Variables ordered by Simulink
            if ~isempty(evalin(mdlWks, [bvarname 'Ploads'])) % If there are loads
                v = evalin(mdlWks, strcat("abs(", bvarname, "lf_op.op_report.Outputs(7).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]) + 1i*", bvarname, "lf_op.op_report.Outputs(8).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]))'"));
                theta = evalin(mdlWks, strcat("angle(", bvarname, "lf_op.op_report.Outputs(7).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]) + 1i*", bvarname, "lf_op.op_report.Outputs(8).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]))'"));
            else
                v = evalin(mdlWks, strcat("abs(", bvarname, "lf_op.op_report.Outputs(5).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]) + 1i*", bvarname, "lf_op.op_report.Outputs(6).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]))'"));
                theta = evalin(mdlWks, strcat("angle(", bvarname, "lf_op.op_report.Outputs(5).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]) + 1i*", bvarname, "lf_op.op_report.Outputs(6).y([", num2str(evalin(mdlWks, ...
                    strcat(bvarname, Generators{i}, "s_buses'"))), "]))'"));
            end
            P = evalin(mdlWks, strcat(bvarname, "lf_op.op_report.Outputs(1).y([", num2str(evalin(mdlWks, ...
                strcat(bvarname, Generators{i}, "_idx'"))), "])'"));
            Q = evalin(mdlWks, strcat(bvarname, "lf_op.op_report.Outputs(2).y([", num2str(evalin(mdlWks, ...
                strcat(bvarname, Generators{i}, "_idx'"))), "])'"));
            if exist(['fun_init_' Generators{i}], 'file')
                feval(['fun_init_' Generators{i}], [block '/Dynamic/Generators/' Generators{i}], v', theta', P', Q', block)
            else
                waitbar(0, wb, ['Initialising ' Generators{i} 's']);
                h = new_system; % Create new simulink model (will not be saved)
                tempWks = get_param(h, 'ModelWorkspace'); % Get temporal model workspace
                assignin(tempWks, 'bvarname', bvarname) % Assign to temporal model workspace
                assignin(tempWks, [bvarname 'lf_op'], evalin(mdlWks, [bvarname 'lf_op'])) % Assign to temporal model workspace
                assignin(tempWks, [Generators{i} 's_buses'], evalin(mdlWks, [bvarname Generators{i} 's_buses']))
                assignin(tempWks, 'buses_lf', evalin(mdlWks, [bvarname 'buses_lf']))
                assignin(tempWks, [Generators{i} '_idx'], evalin(mdlWks, [bvarname Generators{i} '_idx']))
                try
                assignin(tempWks, 'fref', evalin(mdlWks, [bvarname 'fref']))
                catch ME
                end
                assignin(tempWks, 'Nbuses', evalin(mdlWks, [bvarname 'Nbuses']))
                assignin(tempWks, 'gen_static_params', evalin(mdlWks, [bvarname 'gen_static_params']))
                assignin(tempWks, 'Sb', evalin('base', get_param(block, 'Sb')))
                assignin(tempWks, 'iPs_lf', evalin(mdlWks, [bvarname 'iPs_lf']))
                assignin(tempWks, 'iQs_lf', evalin(mdlWks, [bvarname 'iQs_lf']))
                assignin(tempWks, 'Ps_lf', evalin(mdlWks, [bvarname 'Ps_lf']))
                assignin(tempWks, 'Qs_lf', evalin(mdlWks, [bvarname 'Qs_lf']))
                assignin(tempWks, 'omega_base', omega_base)
                assignin(tempWks, 'From_buses', evalin(mdlWks, [bvarname 'From_buses']))
                subsystem_name = get_param(h, 'Name'); % Get name of new model
                add_block('simulink/Signal Routing/Goto', [subsystem_name '/omega_ref'], ...
                    'GotoTag', 'omega_ref', 'TagVisibility', 'global'); % Add omega_ref goto connected to 1 pu
                add_block('simulink/Signal Routing/Goto', [subsystem_name '/omega_base_goto'], ...
                    'GotoTag', 'omega_base', 'TagVisibility', 'global'); % Add omega_base goto
                add_block('simulink/Signal Routing/Selector', [subsystem_name '/selector1'], ...
                    'NumberOfDimensions', '2', 'IndexOptionArray', {'Index vector (dialog)', 'Select all'}, ...
                    'IndexParamArray', {'[1:4]','1'});
                add_block('simulink/Signal Routing/Selector', [subsystem_name '/selector2'], ...
                    'NumberOfDimensions', '2', 'IndexOptionArray', {'Index vector (dialog)', 'Select all'}, ...
                    'IndexParamArray', {'6','1'});
                add_block('simulink/Matrix Operations/Matrix Concatenate', [subsystem_name '/Matrixconcatenate'], ...
                    'NumInputs', '3', 'ConcatenateDimension', '1');
                add_block('simulink/Sources/Constant', [subsystem_name '/one'], ...
                    'Value', 'ones(1,length(omega_base))', 'VectorParams1D', 'off');
                add_block('simulink/Sources/Constant', [subsystem_name '/omega_base'], ...
                    'Value', 'omega_base', 'VectorParams1D', 'off');
                add_block('simulink/Sources/Constant', [subsystem_name '/onev'], ...
                    'Value', 'ones(1,length(buses_lf))', 'VectorParams1D', 'off');
                add_line(subsystem_name, 'one/1', 'omega_ref/1')
                add_line(subsystem_name, 'omega_base/1', 'omega_base_goto/1')
                add_block([block '/Dynamic/Generators/' Generators{i}], [subsystem_name '/' Generators{i}], 'LinkStatus', 'none')
                bl = find_system([subsystem_name '/' Generators{i}], ...
                    'LookUnderMasks', 'on' ,'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
                for ii = 1:length(bl)
                    set_param(bl{ii}, 'InitialCondition', ['ones(1, length(' Generators{i} 's_buses))']) % Initial variables to 1
                end
                waitbar(0.8,wb);
                % Add inputs (voltage of bus of connection)
                in_v = [real(v.*exp(1i*theta))'; imag(v.*exp(1i*theta))'; zeros(1, length(theta))];
                assignin(tempWks, 'in_v', in_v)
                add_block('simulink/Sources/Constant', [subsystem_name '/v'], 'Value', 'in_v');
                elec_in = ones(6, length(P));
                elec_in(1,:) = 1; % fsp
                elec_in(2,:) = P; % Psp
                elec_in(3,:) = v; % vsp
                elec_in(4,:) = Q; % Qsp
                assignin(tempWks, 'elec_in', elec_in)
                outports_block = get_param([subsystem_name '/' Generators{i}],'Ports');
                outports_block = outports_block(2);
    
                add_block('simulink/Sources/Constant', [subsystem_name '/elec'], 'Value', 'elec_in');
                add_line(subsystem_name, 'v/1', [Generators{i} '/2'])
    	        add_block('simulink/Sinks/Out1', [subsystem_name '/' Generators{i} '/out_trim1'])
                add_line([subsystem_name '/' Generators{i}], 'trim_output/1', 'out_trim1/1')
                add_block('simulink/Sinks/Out1', [subsystem_name '/' Generators{i} '/out_trim2'])
                add_line([subsystem_name '/' Generators{i}], 'trim_output/2', 'out_trim2/1')
                add_block('simulink/Sinks/Out1', [subsystem_name '/out_trim1'])
                add_line(subsystem_name, [Generators{i} '/' num2str(outports_block+1)], 'out_trim1/1')
                add_block('simulink/Sinks/Out1', [subsystem_name '/out_trim2'])
                add_line(subsystem_name, [Generators{i} '/' num2str(outports_block+2)], 'out_trim2/1')
    
                % Connect elec variables bus
                add_line(subsystem_name, 'elec/1', 'selector1/1')
                add_line(subsystem_name, 'elec/1', 'selector2/1')
                add_line(subsystem_name, 'selector1/1', 'Matrixconcatenate/1')
                add_line(subsystem_name, 'onev/1', 'Matrixconcatenate/2')
                add_line(subsystem_name, 'selector2/1', 'Matrixconcatenate/3')
                add_line(subsystem_name, 'Matrixconcatenate/1', [Generators{i} '/1'])
               
                if length(theta) > 10 % If more than N units, initialise the first one, then use that point as initial for the others
    
                    % The first one
                    assignin(tempWks, [Generators{i} 's_buses'], evalin(mdlWks, [bvarname Generators{i} 's_buses(1)']))
                    assignin(tempWks, 'buses_lf', '1')
                    assignin(tempWks, [Generators{i} '_idx'], evalin(mdlWks, [bvarname Generators{i} '_idx(1)']))
                    assignin(tempWks, [Generators{i} 's_toM'], evalin(mdlWks, [bvarname Generators{i} 's_toM(1)']))
    
                    % Add inputs (voltage of bus of connection)
                    in_v = [real(v(1).*exp(1i*theta(1)))'; imag(v(1).*exp(1i*theta(1)))'; zeros(1, length(theta(1)))];
                    assignin(tempWks, 'in_v', in_v)
                    elec_in = ones(6, length(P(1)));
                    elec_in(1,:) = 1; % fsp
                    elec_in(2,:) = P(1); % Psp
                    elec_in(3,:) = v(1); % vsp
                    elec_in(4,:) = Q(1); % Qsp
                    assignin(tempWks, 'elec_in', elec_in)
        
                    opspec = operspec(subsystem_name); % Create opspec
                    vd = real(v(1).*exp(1i.*theta(1)));
                    vq = imag(v(1).*exp(1i.*theta(1)));
                    opspec.Outputs(1).Known = true; % iod
                    opspec.Outputs(1).y = (vd.*P(1)+vq.*Q(1))./(v(1).^2);
                    opspec.Outputs(2).Known = true; % ioq
                    opspec.Outputs(2).y = (vq.*P(1)-vd.*Q(1))./(v(1).^2);
                    % optionop = findopOptions('OptimizerType', 'graddescent-elim', 'OptimizationOptions', optimset('PlotFcns', ...
                    %     'optimplotx', "MaxFunEvals", 200*length(theta), "MaxIter", 50));
                    [op.op, op.op_report] = findop(subsystem_name, opspec); % Find operating point
    
                    bl = find_system([subsystem_name '/' Generators{i}], ...
                        'LookUnderMasks', 'on' ,'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
                    for ii = 1:length(bl)
                        set_param(bl{ii}, 'InitialCondition', [num2str(op.op.States(strcmp({op.op.States.Block}, bl{ii})).x) '*ones(1, length(' Generators{i} 's_buses))']) % Initial variables to previous
                    end
                end
    
                % With that initial point, findop for all
    
                assignin(tempWks, [Generators{i} 's_buses'], evalin(mdlWks, [bvarname Generators{i} 's_buses']))
                buses_lf = evalin(mdlWks, [bvarname 'buses_lf']);
                assignin(tempWks, 'buses_lf', buses_lf)
                block_idx = evalin(mdlWks, [bvarname Generators{i} '_idx']);
                assignin(tempWks, [Generators{i} '_idx'], block_idx)
                block_toM = evalin(mdlWks, [bvarname Generators{i} 's_toM']);
                assignin(tempWks, [Generators{i} 's_toM'], block_toM)
                Nbuses = evalin(mdlWks, [bvarname 'Nbuses']);
                assignin(tempWks, 'Nbuses', Nbuses)
                assignin(tempWks, 'Ps_lf', evalin(mdlWks, [bvarname 'Ps_lf']))
                assignin(tempWks, 'Qs_lf', evalin(mdlWks, [bvarname 'Qs_lf']))
                assignin(tempWks, 'From_buses', evalin(mdlWks, [bvarname 'From_buses']))
        
                in_v = [real(v.*exp(1i*theta))'; imag(v.*exp(1i*theta))'; zeros(1, length(theta))];
                assignin(tempWks, 'in_v', in_v)
                elec_in = ones(6, length(buses_lf));
                elec_in(1, block_idx) = 1; % fsp
                elec_in(2, block_idx) = P; % Psp
                elec_in(3, block_idx) = v; % vsp
                elec_in(4, block_idx) = Q; % Qsp
                assignin(tempWks, 'elec_in', elec_in)
    
               
                opspec = operspec(subsystem_name); % Create opspec
                vd = real(v.*exp(1i.*theta));
                vq = imag(v.*exp(1i.*theta));
                opspec.Outputs(1).Known = true(1,  length(theta)); % iod
                opspec.Outputs(1).y = (vd.*P+vq.*Q)./(v.^2);
                opspec.Outputs(2).Known = true(1,  length(theta)); % ioq
                opspec.Outputs(2).y = (vq.*P-vd.*Q)./(v.^2);
                optionop = findopOptions('OptimizerType', 'graddescent-elim', 'OptimizationOptions', optimset("MaxFunEvals", 2000*length(theta), "MaxIter", 50,'TolFun',1e-6,'TolX',1e-6,'PlotFcns',@optimplotx));
                optionop.OptimizationOptions.Algorithm = 'interior-point';
                tic
                [op.op, op.op_report] = findop(subsystem_name, opspec, optionop); % Find operating point
                toc

                if ~strcmp(op.op_report.TerminationString, 'Operating point specifications were successfully met.')
                    assignin(mdlWks, [bvarname 'init_OK_bool'], false)
                end
                bl = find_system([subsystem_name '/' Generators{i}], ...
                    'LookUnderMasks', 'on' ,'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
                for ii = 1:length(bl)
                    set_param(bl{ii}, 'InitialCondition', ['ones(1, length(' Generators{i} 's_buses))']) % Initial variables to previous
                end
    
                close % close PlotFcns
                if ~strcmp(op.op_report.TerminationString, "Operating point specifications were successfully met.")
                    errordlg("Generate loadflow report for details", [Generators{i} ' initialisation did not converge'])
                end
                
                assignin(mdlWks, [bvarname 'op_' Generators{i}], op) % Assign to model workspace
                bdclose(subsystem_name) % Close temporal model
                
                % Load the operating point to the model (paths are different)
                for ii = 1:length(op.op.States)
                    newpath = string(strrep({op.op.States(ii).Block}, subsystem_name, [block '/Dynamic/Generators']));
                    set_param(newpath, 'InitialCondition', strcat(bvarname, "op_", Generators{i}, ".op.States(", num2str(ii), ").x'"))
                end
            waitbar(1,wb);
            end
        end
        %% ----------------     DC system    ------------------------------
        [Tbus, ~, ~, TbusDC, ~, Tintconv, ~, ~] = fun_read_mask_tables(gcb);
        if isempty(Tintconv)
            set_param([block '/Dynamic/DCAC converters'], 'Commented', 'on')
            set_param([block '/Dynamic/Grid_DC'], 'Commented', 'on')
        else
            set_param([block '/Dynamic/DCAC converters'], 'Commented', 'off')
            set_param([block '/Dynamic/Grid_DC'], 'Commented', 'off')
            DC_DGs = Tintconv.type;
            DC_intconv = unique(DC_DGs);
            DC = evalin(mdlWks, [bvarname 'DC']);
            bl = find_system([block '/Dynamic/DCAC converters'], ...
                'IncludeCommented', 'on', 'SearchDepth', '2', ...
                'LookUnderMasks', 'on', 'Name', 'trim_output');
            for i = 1:length(bl)
                bl{i} = get_param(bl{i}, 'Parent');
            end
            for i=1:length(bl) 
                if sum(strcmp(DC_DGs, get_param(bl{i}, 'Name'))) == 0
                    set_param(bl{i}, "commented", "on")
                else
                    set_param(bl{i}, "commented", "off")
                end
            end
            DGs = TgenDC.type;
            GeneratorsDC = unique(DGs);

            % Search for depth up to 2 (includes depth 1 and 2)
            all_depth_results = find_system([block '/Dynamic/Generators_DC'], ...
                                'IncludeCommented', 'on', ...
                                'SearchDepth', '2', ...
                                'LookUnderMasks', 'on', ...
                                'Name', 'ref');
            % Search for depth 1 only
            depth_1_results = find_system([block '/Dynamic/Generators_DC'], ...
                              'IncludeCommented', 'on', ...
                              'SearchDepth', '1', ...
                              'LookUnderMasks', 'on', ...
                              'Name', 'ref');

            % Filter out depth 1 results from the all_depth_results
            bl = setdiff(all_depth_results, depth_1_results);
            for i = 1:length(bl)
                bl{i} = get_param(bl{i}, 'Parent');
            end
            for i=1:length(bl) 
                if sum(strcmp(GeneratorsDC, get_param(bl{i}, 'Name'))) == 0
                    set_param(bl{i}, "commented", "on")
                else
                    set_param(bl{i}, "commented", "off")
                end
            end

            if ~isempty(GeneratorsDC) 
                for i=1:length(GeneratorsDC)
                    
                    % Get the buses where DC units are connected
                    DGs_typei_location = find(strcmp(DGs, GeneratorsDC{i}));
                    aux = zeros(1, length(DGs_typei_location));
                    for ii = 1:length(aux)
                        aux(ii) = find(TbusDC.bus_i == TgenDC.bus(DGs_typei_location(ii)));
                    end
                    DC.(GeneratorsDC{i} + "s_buses") = aux;
                    %assignin(mdlWks, [bvarname GeneratorsDC{i} 's_buses'], aux')
        
                    aux_M = zeros(height(TbusDC), length(DGs_typei_location));
                    aux_idx = sub2ind(size(aux_M), aux, 1:length(DGs_typei_location));
                    aux_M(aux_idx) = 1;
                    DC.(GeneratorsDC{i} + "s_toM") = aux_M;
                    DC.(GeneratorsDC{i} + "s_Pg") = TgenDC.Pg(DGs_typei_location(ii));
                    %assignin(mdlWks, [bvarname GeneratorsDC{i} 's_toM'], aux_M)


                    h = new_system; % Create new simulink model (will not be saved)
                    tempWks = get_param(h, 'ModelWorkspace'); % Get temporal model workspace
                    subsystem_name = get_param(h, 'Name'); % Get name of new model
                    add_block([block '/Dynamic/Generators_DC/' GeneratorsDC{i}], [subsystem_name '/' GeneratorsDC{i}], 'LinkStatus', 'none')
                    bl = find_system([subsystem_name '/' GeneratorsDC{i}], ...
                        'LookUnderMasks', 'on' ,'MatchFilter', @Simulink.match.activeVariants, 'BlockType', 'Integrator');
                    for ii = 1:length(bl)
                        set_param(bl{ii}, 'InitialCondition', ['ones(1, length(DC.' GeneratorsDC{i} 's_buses))']) % Initial variables to 1
                    end
                    % Add inputs (voltage of bus of connection)
                    
                    Sb = evalin('base', get_param(block, 'Sb'));
                    assignin(tempWks, 'Sb', Sb)
                    add_block('simulink/Sources/Constant', [subsystem_name '/v'], 'Value', 'DC.v');
                    add_line(subsystem_name, 'v/1', [GeneratorsDC{i} '/1'])
                    add_block('simulink/Sinks/Out1', [subsystem_name '/out_trim1'])
                    add_line(subsystem_name, [GeneratorsDC{i} '/1'], 'out_trim1/1')
                    DCi.v = DC.v(DC.(GeneratorsDC{i} + "s_buses"));
                    DCi.Nbuses = length(DCi.v);
                    DCi.(GeneratorsDC{i} + "s_Pg") = DC.(GeneratorsDC{i} + "s_Pg");
                    DCi.(GeneratorsDC{i} + "s_buses") = DC.(GeneratorsDC{i} + "s_buses");
                    assignin(tempWks, 'DC', DCi)
                    opspec = operspec(subsystem_name); % Create opspec
                    opspec.Outputs(1).Known = true(1,  length(DCi.v)); % i
                    opspec.Outputs(1).y = DC.(GeneratorsDC{i} + "s_Pg")/Sb./DCi.v;
                    optionop = findopOptions('OptimizerType', 'graddescent-elim', 'OptimizationOptions', ...
                        optimset("MaxFunEvals", 200*length(theta), "MaxIter", 50));
                    clear op
                    [op.op, op.op_report] = findop(subsystem_name, opspec, optionop); % Find operating point
                    bdclose(subsystem_name) % Close temporal model
                    assignin(mdlWks, [bvarname 'DCop_' GeneratorsDC{i}], op)
                    % Load the operating point to the model (paths are different)
                    for ii = 1:length(op.op.States)
                        newpath = string(strrep({op.op.States(ii).Block}, subsystem_name, [block '/Dynamic/Generators_DC']));
                        set_param(newpath, 'InitialCondition', strcat(bvarname, "DCop_", GeneratorsDC{i}, ".op.States(", num2str(ii), ").x'"))
                    end
                end
                assignin(mdlWks, [bvarname 'DC'], DC)
            end
            
            for i=1:length(DC_intconv)
                % Get the buses where units are connected
                DGs_typei_location = find(strcmp(DC_DGs, DC_intconv{i}));
                aux = zeros(length(DGs_typei_location), 1);
                for ii = 1:length(aux)
                    aux(ii) = find(Tbus.bus_i == Tintconv.tbus(DGs_typei_location(ii)));
                end
                aux_M = zeros(height(Tbus), length(DGs_typei_location));
                aux_idx = sub2ind(size(aux_M), aux', 1:length(DGs_typei_location));
                aux_M(aux_idx) = 1;

                % Get the buses where units are connected in DC
                DGs_typei_location = find(strcmp(DC_DGs, DC_intconv{i}));
                aux_DC = zeros(length(DGs_typei_location), 1);
                for ii = 1:length(aux_DC)
                    aux_DC(ii) = find(TbusDC.bus_i == Tintconv.fbus(DGs_typei_location(ii)));
                end
                aux_DC_M = zeros(height(TbusDC), length(DGs_typei_location));
                aux_DC_idx = sub2ind(size(aux_DC_M), aux_DC', 1:length(DGs_typei_location));
                aux_DC_M(aux_DC_idx) = 1;
                DC.(DC_intconv{i} + "s_toM") = aux_M;
                DC.(DC_intconv{i} + "s_DC_toM") = aux_DC_M;
                DC.(DC_intconv{i} + "s_buses") = aux;
                DC.(DC_intconv{i} + "s_DC_buses") = aux_DC;
                dparams.ratio = Tintconv.ratio(DGs_typei_location)';
                dparams.ratio(dparams.ratio == 0) = 1;
                dparams.type = Tintconv.CONV(DGs_typei_location)';
                dparams.Rf = Tintconv.r(DGs_typei_location)';
                dparams.Lf = Tintconv.x(DGs_typei_location)';
                dparams.alpha1 = Tintconv.ALPHA1(DGs_typei_location)';
                dparams.alpha2 = Tintconv.ALPHA2(DGs_typei_location)';
                dparams.alpha3 = Tintconv.ALPHA3(DGs_typei_location)';
                dparams.Gsw = Tintconv.GSW(DGs_typei_location)';
                DC.(DC_intconv{i} + "_dparams") = dparams;
                assignin(mdlWks, [bvarname 'DC'], DC)
                % Initialise
                intconv_PF = evalin(mdlWks, [bvarname 'intconv_PF']);
                v = intconv_PF.vmod(DGs_typei_location);
                theta = deg2rad(intconv_PF.vangle(DGs_typei_location));
                P = -intconv_PF.P(DGs_typei_location);
                Q = -intconv_PF.Q(DGs_typei_location);
                feval(['fun_init_intconv_' DC_intconv{i}], [block '/Dynamic/DCAC converters/' DC_intconv{i}], v, theta, P, Q, block)
            end
        end
        close(wb)
        catch ME
            if length(ME.cause) > 1
                message_multiple_causes = strings(length(ME.cause), 1);
                for i=1:length(ME.cause)
                    message_multiple_causes(i) = ME.cause{i}.message;
                end
                str_name = string({ME.stack.name}');
                str_line = string({ME.stack.line}');
                errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
                errordlg(char(errortext));
                try
                    fig = uifigure('Name', 'Error', 'Position', [100 100 400 50*length(ME.cause)]);
                    uilabel(fig, 'Text', message_multiple_causes, 'Interpreter', 'html', 'Position', [10 50 380 50*length(ME.cause)]);
                    uibutton(fig, 'Text', 'OK', 'Position', [170 10 60 30], 'ButtonPushedFcn', @(btn,event) close(fig));
                catch
                end
            end
            str_name = string({ME.stack.name}');
            str_line = string({ME.stack.line}');
            errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
            errordlg(char(errortext));
        end
    end
end