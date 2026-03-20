function fun_lf_report_callback(gcb)
try
    block = gcb;
    [Tbus, Tgen, Tline, TbusDC, TlineDC, Tintconv, ...
        Tbuscomplete, Tlinecomplete, TFb, TgenDC, Tgencomplete] = fun_read_mask_tables(block);
    if ~isempty(Tintconv)
        Tlinecomplete_withintconv = [Tlinecomplete(:, 1:13); Tintconv(:, 2:14)];
    else
        Tlinecomplete_withintconv = Tlinecomplete;
    end
    mdlWks = get_param(bdroot(gcb), 'ModelWorkspace');
    bvarname = get_param(gcb, 'bvarname');
    PF_branch_complete = evalin(mdlWks, [bvarname 'PF_branch_complete']);
    if size(PF_branch_complete, 2) == 37
        linePF = PF_branch_complete(:,14);
        lineQF = PF_branch_complete(:,15);
        linePT = PF_branch_complete(:,16);
        lineQT = PF_branch_complete(:,17);
    end
    Ggrid = graph(string(Tlinecomplete_withintconv.fbus), string(Tlinecomplete_withintconv.tbus), abs(Tlinecomplete_withintconv.x + 1i*Tlinecomplete_withintconv.x));
    [~, idx_nodes] = ismember(str2double(Ggrid.Nodes.Name), Tbuscomplete.bus_i);
    Ggrid.Nodes.Pload = Tbuscomplete.Pd(idx_nodes);
    Ggrid.Nodes.Qload = Tbuscomplete.Qd(idx_nodes);
    Pgen_bus = zeros(height(Tbuscomplete), 1);
    Qgen_bus = zeros(height(Tbuscomplete), 1);
    for i = 1:height(Tbuscomplete)
        Pgen_bus(i) = sum(Tgencomplete.Pg(Tgencomplete.bus == Tbuscomplete.bus_i(i)));
        Qgen_bus(i) = sum(Tgencomplete.Qg(Tgencomplete.bus == Tbuscomplete.bus_i(i)));
    end
    for i = 1:height(Ggrid.Nodes)
        Ggrid.Nodes.Pgen(i) = Pgen_bus(idx_nodes(i));
        Ggrid.Nodes.Qgen(i) = Qgen_bus(idx_nodes(i));
    end
    if ~isempty(TFb)
        areas = unique(TFb.Area);
    else
        areas = [];
    end
    figure(Position=[100 100 600 500])
    if any(Tbus.Latitude ~= 0)
        [~, from] = ismember(Tlinecomplete.fbus, Tbuscomplete.bus_i);
        [~, to] = ismember(Tlinecomplete.tbus, Tbuscomplete.bus_i);

        % --- 2. Compute Midpoints for Arrows ---
        mid_lat = (Tbuscomplete.Latitude(from) + Tbuscomplete.Latitude(to)) / 2;
        mid_lon = (Tbuscomplete.Longitude(from) + Tbuscomplete.Longitude(to)) / 2;

        % --- 3. Plot the Cables (Lines) with Flow Labels ---
        for i = 1:height(Tlinecomplete)
            % Coordinates (Include midpoint explicitly to fix arrowhead detachment on zoom)
            lats = [Tbuscomplete.Latitude(from(i)), mid_lat(i), Tbuscomplete.Latitude(to(i))];
            lons = [Tbuscomplete.Longitude(from(i)), mid_lon(i), Tbuscomplete.Longitude(to(i))];

            % Draw the line
            geoplot(lats, lons, 'white-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            hold on
        end

        % --- 4. Compute Power Flow Direction and Branch Losses ---
        P_sending = zeros(height(Tlinecomplete), 1);
        P_receiving = zeros(height(Tlinecomplete), 1);
        Ploss = zeros(height(Tlinecomplete), 1);
        Ploss_pct = zeros(height(Tlinecomplete), 1);
        dlat = zeros(height(Tlinecomplete), 1);
        dlon = zeros(height(Tlinecomplete), 1);
        from_bus_flow = zeros(height(Tlinecomplete), 1);
        to_bus_flow = zeros(height(Tlinecomplete), 1);

        for i = 1:height(Tlinecomplete)
            if linePF(i) >= 0
                P_sending(i) = linePF(i);
                P_receiving(i) = -linePT(i);
                dlat(i) = Tbuscomplete.Latitude(to(i)) - Tbuscomplete.Latitude(from(i));
                dlon(i) = Tbuscomplete.Longitude(to(i)) - Tbuscomplete.Longitude(from(i));
                from_bus_flow(i) = Tlinecomplete.fbus(i);
                to_bus_flow(i) = Tlinecomplete.tbus(i);
            else
                P_sending(i) = linePT(i);
                P_receiving(i) = -linePF(i);
                dlat(i) = Tbuscomplete.Latitude(from(i)) - Tbuscomplete.Latitude(to(i));
                dlon(i) = Tbuscomplete.Longitude(from(i)) - Tbuscomplete.Longitude(to(i));
                from_bus_flow(i) = Tlinecomplete.tbus(i);
                to_bus_flow(i) = Tlinecomplete.fbus(i);
            end
            Ploss(i) = linePF(i) + linePT(i);

            if abs(P_sending(i)) > 0
                Ploss_pct(i) = Ploss(i) / abs(P_sending(i)) * 100;
            end
        end
        angles = atan2d(dlat, dlon);

        idx_east = angles >= -45 & angles < 45;
        idx_north = angles >= 45 & angles < 135;
        idx_west = angles >= 135 | angles < -135;
        idx_south = angles >= -135 & angles < -45;

        % --- 5. Plot Branches Flow Arrows ---
        markers = {'>', '^', '<', 'v'};
        inds = {idx_east, idx_north, idx_west, idx_south};
        for j = 1:4
            idx_sub = inds{j};
            if any(idx_sub)
                s_arr = geoscatter(mid_lat(idx_sub), mid_lon(idx_sub), 60, markers{j}, 'filled', ...
                    'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'k');
                s_arr.DataTipTemplate.DataTipRows = [ ...
                    dataTipTextRow('From Bus', from_bus_flow(idx_sub)), ...
                    dataTipTextRow('To Bus', to_bus_flow(idx_sub)), ...
                    dataTipTextRow('Sending P (MW)', P_sending(idx_sub)), ...
                    dataTipTextRow('Receiving P (MW)', P_receiving(idx_sub)), ...
                    dataTipTextRow('P Loss', string(round(Ploss(idx_sub),4)) + " MW (" + string(round(Ploss_pct(idx_sub),2)) + "%)") ...
                    ];
            end
        end

        % --- 5. Plot the Nodes (Buses) with Custom Datatips ---
        % We use geoscatter because it supports 'DataTipTemplate'
        s = geoscatter(Tbuscomplete.Latitude, Tbuscomplete.Longitude, 50, 'o', 'filled', ...
            'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k');

        % Configure the Datatip for buses
        % Configure the Datatip for buses
        s.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('Bus ID', Tbuscomplete.bus_i), ...
            dataTipTextRow('Lat', Tbuscomplete.Latitude), ...
            dataTipTextRow('Lon', Tbuscomplete.Longitude), ...
            dataTipTextRow('P gen (MW)', Pgen_bus), ...
            dataTipTextRow('Q gen (MVAr)', Qgen_bus), ...
            dataTipTextRow('P load (MW)', Tbuscomplete.Pd), ...
            dataTipTextRow('Q load (MVAr)', Tbuscomplete.Qd), ...
            dataTipTextRow('Voltage (pu)', Tbuscomplete.Vm), ...
            dataTipTextRow('Base kV', Tbuscomplete.baseKV) ...
            ];

    else
        h = plot(Ggrid, 'Layout', 'force','UseGravity',true,'DisplayName', 'Grid graph');
        highlight(h, (Ggrid.Nodes.Pgen ~= 0 | Ggrid.Nodes.Qgen ~= 0), 'NodeColor', 'red')
        lineDC = [startsWith(Tlinecomplete.Complexity, 'DC '); false(height(Tintconv), 1)];
        lineDCAC = [false(height(Tlinecomplete), 1); true(height(Tintconv), 1)];
        highlight(h, string(Tlinecomplete_withintconv.fbus(lineDC)), string(Tlinecomplete_withintconv.tbus(lineDC)),'LineStyle','--','EdgeColor','green','LineWidth',1.5)
        plotslgnd.p1 = plot(nan, nan, 'green','LineStyle','--','LineWidth',1.5,'DisplayName', 'DC grid');
        highlight(h, string(Tlinecomplete_withintconv.fbus(lineDCAC)), string(Tlinecomplete_withintconv.tbus(lineDCAC)),'LineStyle','-.','EdgeColor','black','LineWidth',1.5)
        plotslgnd.p2 = plot(nan, nan, 'black','LineStyle','-.','LineWidth',1.5,'DisplayName', 'DC/AC converters');
        colors_area = ["magenta" "white"];
        lines_style_area = ["-" "-"];
        for i = 1:length(areas)
            areas(i);
            larea = ismember(Tlinecomplete_withintconv.fbus, Tbuscomplete.bus_i(Tbuscomplete.area == areas(i))) & ...
                ismember(Tlinecomplete_withintconv.tbus, Tbuscomplete.bus_i(Tbuscomplete.area == areas(i)));
            color_area = colors_area(mod(i, length(colors_area)) + 1);
            line_style_area = lines_style_area(mod(i, length(lines_style_area)) + 1);
            highlight(h, string(Tlinecomplete_withintconv.fbus(larea)), string(Tlinecomplete_withintconv.tbus(larea)),'EdgeColor',color_area,'LineWidth',1.5, 'LineStyle', line_style_area)
            plotslgnd.("p" + 2 + i) = plot(nan, nan, color_area,'LineWidth',1.5,'LineStyle', line_style_area,'DisplayName', ['AC area ' num2str(areas(i))]);
        end
        legend('Location','best');
        hold on
        h.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('Bus ID', Ggrid.Nodes.Name), ...
            dataTipTextRow('Load P (MW)', Ggrid.Nodes.Pload), ...
            dataTipTextRow('Load Q (MVAr)', Ggrid.Nodes.Qload), ...
            dataTipTextRow('Gen P (MW)', Ggrid.Nodes.Pgen), ...
            dataTipTextRow('Gen Q (MVAr)', Ggrid.Nodes.Qgen), ...
            dataTipTextRow('V (pu)', Tbuscomplete.Vm) ...
            ];
    end


catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end