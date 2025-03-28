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
    Ggrid = graph(string(Tlinecomplete_withintconv.fbus), string(Tlinecomplete_withintconv.tbus), abs(Tlinecomplete_withintconv.x + 1i*Tlinecomplete_withintconv.x));
    [~, idx_nodes] = ismember(str2double(Ggrid.Nodes.Name), Tbuscomplete.bus_i);
    Ggrid.Nodes.Pload = Tbuscomplete.Pd(idx_nodes);
    Ggrid.Nodes.Qload = Tbuscomplete.Qd(idx_nodes);
    for i = 1:height(Ggrid.Nodes)
        Ggrid.Nodes.Pgen(i) = sum(Tgencomplete.Pg(Tgencomplete.bus == Tbuscomplete.bus_i(idx_nodes(i))));
        Ggrid.Nodes.Qgen(i) = sum(Tgencomplete.Qg(Tgencomplete.bus == Tbuscomplete.bus_i(idx_nodes(i))));
    end
    if ~isempty(TFb)
        areas = unique(TFb.Area);
    else
        areas = [];
    end
    figure(Position=[100 100 600 500])
    h = plot(Ggrid, 'Layout', 'force','UseGravity',true,'DisplayName', 'Grid graph');
    hold on
    h.DataTipTemplate.DataTipRows(2) = [];
    h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Load P (MW)', Ggrid.Nodes.Pload);
    h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Load Q (MVAr)', Ggrid.Nodes.Qload);
    h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Gen P (MW)', Ggrid.Nodes.Pgen);
    h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Gen Q (MVAr)', Ggrid.Nodes.Qgen);
    h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('V (pu)', Tbuscomplete.Vm);
    highlight(h, (Ggrid.Nodes.Pgen ~= 0 | Ggrid.Nodes.Qgen ~= 0), 'NodeColor', 'red')
    lineDC = [startsWith(Tlinecomplete.Complexity, 'DC '); false(height(Tintconv), 1)];
    lineDCAC = [false(height(Tlinecomplete), 1); true(height(Tintconv), 1)];
    highlight(h, string(Tlinecomplete_withintconv.fbus(lineDC)), string(Tlinecomplete_withintconv.tbus(lineDC)),'LineStyle','--','EdgeColor','green','LineWidth',1.5)
    plotslgnd.p1 = plot(nan, nan, 'green','LineStyle','--','LineWidth',1.5,'DisplayName', 'DC grid');
    highlight(h, string(Tlinecomplete_withintconv.fbus(lineDCAC)), string(Tlinecomplete_withintconv.tbus(lineDCAC)),'LineStyle','-.','EdgeColor','black','LineWidth',1.5)
    plotslgnd.p2 = plot(nan, nan, 'black','LineStyle','-.','LineWidth',1.5,'DisplayName', 'DC/AC converters');
    colors_area = ["magenta" "blue"];
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
catch ME
    str_name = string({ME.stack.name}');
    str_line = string({ME.stack.line}');
    errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
    errordlg(char(errortext));
end
end