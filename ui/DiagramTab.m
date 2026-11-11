classdef DiagramTab < handle
    %DIAGRAMTAB Interactive AC/DC network topology, flows & device inspector tab.

    properties
        app
        parent
        ddLayout
        ddColorBy
        chkFlows
        ddBusInspect
        axNet
        txtInspect
        nodeCoords double = zeros(0, 2)
    end

    methods
        function obj = DiagramTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [2, 1]);
            gl.RowHeight = {40, '1x'};
            gl.Padding = [10, 10, 10, 10];
            gl.RowSpacing = 8;

            % Toolbar
            bar = uigridlayout(gl, [1, 8]);
            bar.ColumnWidth = {95, 55, 115, 72, 165, 130, 155, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uibutton(bar, 'Text', 'Refresh', ...
                'ButtonPushedFcn', @(~,~) obj.refresh());
            uilabel(bar, 'Text', 'Layout:', 'HorizontalAlignment', 'right');
            obj.ddLayout = uidropdown(bar, ...
                'Items', {'force', 'layered', 'circle', 'coordinates'}, ...
                'Value', 'force', ...
                'ValueChangedFcn', @(~,~) obj.refresh());
            uilabel(bar, 'Text', 'Color by:', 'HorizontalAlignment', 'right');
            obj.ddColorBy = uidropdown(bar, ...
                'Items', {'Domain (AC / DC)', 'Voltage Vm (pu)', 'Angle Va (deg)'}, ...
                'Value', 'Domain (AC / DC)', ...
                'ValueChangedFcn', @(~,~) obj.refresh());
            obj.chkFlows = uicheckbox(bar, ...
                'Text', 'Show flow labels', ...
                'Value', true, ...
                'ValueChangedFcn', @(~,~) obj.refresh());
            uibutton(bar, 'Text', 'Device Parameters...', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'ButtonPushedFcn', @(~,~) obj.openDeviceForSelectedBus());
            uilabel(bar, 'Text', '');

            % Split canvas + right-side inspector
            body = uigridlayout(gl, [1, 2]);
            body.ColumnWidth = {'1x', 320};
            body.Padding = [0, 0, 0, 0];
            body.ColumnSpacing = 10;

            obj.axNet = uiaxes(body);
            obj.axNet.Color = [1, 1, 1];
            obj.axNet.Box = 'on';

            pInsp = uipanel(body, 'Title', 'Bus & Element Inspector');
            gInsp = uigridlayout(pInsp, [2, 1]);
            gInsp.RowHeight = {32, '1x'};
            gInsp.Padding = [8, 8, 8, 8];

            selBar = uigridlayout(gInsp, [1, 2]);
            selBar.ColumnWidth = {80, '1x'};
            selBar.Padding = [0, 0, 0, 0];
            uilabel(selBar, 'Text', 'Inspect bus:');
            obj.ddBusInspect = uidropdown(selBar, ...
                'Items', {'-'}, ...
                'Value', '-', ...
                'ValueChangedFcn', @(~,~) obj.updateInspector());

            obj.txtInspect = uitextarea(gInsp, ...
                'Editable', 'off', ...
                'FontName', 'Consolas', ...
                'FontSize', 11);
        end

        function refresh(obj)
            sess = obj.app.session;
            TBus = sess.tables.Bus;
            TLine = sess.tables.Line;
            TGen = sess.tables.Gen;
            TIC = sess.tables.IntConv;

            cla(obj.axNet);
            if isempty(TBus) || height(TBus) == 0
                title(obj.axNet, 'No buses loaded');
                obj.ddBusInspect.Items = {'-'};
                obj.ddBusInspect.Value = '-';
                obj.txtInspect.Value = {'(Load a network case to inspect buses and branches)'};
                return;
            end

            nB = height(TBus);
            busIds = round(TBus.bus_i(:));
            busItems = arrayfun(@(b) sprintf('Bus %d', b), busIds, 'UniformOutput', false);
            prevBus = char(obj.ddBusInspect.Value);
            obj.ddBusInspect.Items = busItems;
            if ismember(prevBus, busItems)
                obj.ddBusInspect.Value = prevBus;
            else
                obj.ddBusInspect.Value = busItems{1};
            end

            % Compute 2D coordinates using MATLAB graph layout or geographical coords
            coords = obj.computeLayout(TBus, TLine, TIC);
            obj.nodeCoords = coords;

            hold(obj.axNet, 'on');
            isDCBus = startsWith(strtrim(TBus.Complexity_bus), "DC", 'IgnoreCase', true);

            % 1. Draw AC & DC lines
            if ~isempty(TLine) && height(TLine) > 0
                for k = 1:height(TLine)
                    if TLine.status(k) <= 0, continue; end
                    i1 = find(busIds == round(TLine.fbus(k)), 1);
                    i2 = find(busIds == round(TLine.tbus(k)), 1);
                    if isempty(i1) || isempty(i2), continue; end
                    isDCL = startsWith(strtrim(TLine.Complexity(k)), "DC", 'IgnoreCase', true);
                    if isDCL
                        plot(obj.axNet, [coords(i1,1), coords(i2,1)], [coords(i1,2), coords(i2,2)], ...
                            '--', 'LineWidth', 2.0, 'Color', [0.85, 0.42, 0.10]);
                    else
                        plot(obj.axNet, [coords(i1,1), coords(i2,1)], [coords(i1,2), coords(i2,2)], ...
                            '-', 'LineWidth', 1.9, 'Color', [0.22, 0.44, 0.72]);
                    end
                end
            end

            % 2. Draw Interlinking Converters
            if ~isempty(TIC) && height(TIC) > 0
                for k = 1:height(TIC)
                    if TIC.status(k) <= 0, continue; end
                    i1 = find(busIds == round(TIC.fbus(k)), 1);
                    i2 = find(busIds == round(TIC.tbus(k)), 1);
                    if isempty(i1) || isempty(i2), continue; end
                    plot(obj.axNet, [coords(i1,1), coords(i2,1)], [coords(i1,2), coords(i2,2)], ...
                        '-.', 'LineWidth', 2.6, 'Color', [0.58, 0.22, 0.72]);
                    xm = 0.5 * (coords(i1,1) + coords(i2,1));
                    ym = 0.5 * (coords(i1,2) + coords(i2,2));
                    text(obj.axNet, xm, ym + 0.04, char(TIC.type(k)), ...
                        'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.45, 0.14, 0.58], ...
                        'HorizontalAlignment', 'center', 'Interpreter', 'none');
                end
            end

            % 3. Optional solved active power flow labels
            if obj.chkFlows.Value && ~isempty(sess.pfResult) && isfield(sess.pfResult, 'branchFlows')
                bf = sess.pfResult.branchFlows;
                for k = 1:height(bf)
                    i1 = find(busIds == round(bf.fbus(k)), 1);
                    i2 = find(busIds == round(bf.tbus(k)), 1);
                    if isempty(i1) || isempty(i2), continue; end
                    xm = 0.5 * (coords(i1,1) + coords(i2,1));
                    ym = 0.5 * (coords(i1,2) + coords(i2,2));
                    text(obj.axNet, xm, ym - 0.04, sprintf('%.1f MW', bf.Pf_MW(k)), ...
                        'FontSize', 8, 'Color', [0.18, 0.22, 0.28], ...
                        'BackgroundColor', [1, 1, 1, 0.78], 'HorizontalAlignment', 'center');
                end
            end

            % 4. Draw Bus Nodes
            colorMode = char(obj.ddColorBy.Value);
            if strcmp(colorMode, 'Voltage Vm (pu)')
                cVals = TBus.Vm(:);
                acIdx = find(~isDCBus);
                dcIdx = find(isDCBus);
                if ~isempty(acIdx)
                    scatter(obj.axNet, coords(acIdx,1), coords(acIdx,2), 160, cVals(acIdx), ...
                        'filled', 'Marker', 'o', 'MarkerEdgeColor', [0.1, 0.1, 0.1], ...
                        'ButtonDownFcn', @(~, evt) obj.onCanvasClick(evt));
                end
                if ~isempty(dcIdx)
                    scatter(obj.axNet, coords(dcIdx,1), coords(dcIdx,2), 175, cVals(dcIdx), ...
                        'filled', 'Marker', 's', 'MarkerEdgeColor', [0.1, 0.1, 0.1], ...
                        'ButtonDownFcn', @(~, evt) obj.onCanvasClick(evt));
                end
                colormap(obj.axNet, parula);
                colorbar(obj.axNet);
            elseif strcmp(colorMode, 'Angle Va (deg)')
                cVals = TBus.Va(:);
                scatter(obj.axNet, coords(:,1), coords(:,2), 160, cVals, ...
                    'filled', 'MarkerEdgeColor', [0.1, 0.1, 0.1], ...
                    'ButtonDownFcn', @(~, evt) obj.onCanvasClick(evt));
                colormap(obj.axNet, turbo);
                colorbar(obj.axNet);
            else
                colorbar(obj.axNet, 'off');
                acIdx = find(~isDCBus);
                dcIdx = find(isDCBus);
                if ~isempty(acIdx)
                    scatter(obj.axNet, coords(acIdx,1), coords(acIdx,2), 165, ...
                        repmat([0.15, 0.48, 0.85], numel(acIdx), 1), ...
                        'filled', 'Marker', 'o', 'MarkerEdgeColor', [0.08, 0.22, 0.42], ...
                        'ButtonDownFcn', @(~, evt) obj.onCanvasClick(evt));
                end
                if ~isempty(dcIdx)
                    scatter(obj.axNet, coords(dcIdx,1), coords(dcIdx,2), 175, ...
                        repmat([0.88, 0.44, 0.12], numel(dcIdx), 1), ...
                        'filled', 'Marker', 's', 'MarkerEdgeColor', [0.45, 0.20, 0.04], ...
                        'ButtonDownFcn', @(~, evt) obj.onCanvasClick(evt));
                end
            end

            % 5. Node & Generator Badges
            for i = 1:nB
                bId = busIds(i);
                uList = strings(0, 1);
                if ~isempty(TGen) && height(TGen) > 0
                    m = round(TGen.bus) == bId & TGen.status > 0;
                    if any(m)
                        uList = TGen.type(m);
                    end
                end
                if isempty(uList)
                    lbl = sprintf('Bus %d\n(%.3f pu)', bId, TBus.Vm(i));
                else
                    lbl = sprintf('Bus %d [%s]\n(%.3f pu)', bId, strjoin(uList, ','), TBus.Vm(i));
                end
                text(obj.axNet, coords(i,1), coords(i,2) + 0.08, lbl, ...
                    'FontSize', 8.5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                    'Interpreter', 'none');
            end

            obj.axNet.ButtonDownFcn = @(~, evt) obj.onCanvasClick(evt);
            hold(obj.axNet, 'off');
            axis(obj.axNet, 'padded');
            obj.axNet.XTick = [];
            obj.axNet.YTick = [];
            title(obj.axNet, sprintf('%s — %d AC Buses (circles), %d DC Buses (squares)', ...
                sess.caseName, sum(~isDCBus), sum(isDCBus)), 'Interpreter', 'none');

            obj.updateInspector();
        end

        function coords = computeLayout(obj, TBus, TLine, TIC)
            nB = height(TBus);
            mode = char(obj.ddLayout.Value);
            if strcmp(mode, 'coordinates') && any(TBus.Latitude ~= 0 | TBus.Longitude ~= 0)
                x = TBus.Longitude(:);
                y = TBus.Latitude(:);
                rx = max(x) - min(x); if rx == 0, rx = 1; end
                ry = max(y) - min(y); if ry == 0, ry = 1; end
                coords = [(x - min(x)) / rx, (y - min(y)) / ry];
                return;
            end
            if nB == 1
                coords = [0, 0];
                return;
            end

            busIds = round(TBus.bus_i(:));
            sEdges = [];
            tEdges = [];
            if ~isempty(TLine) && height(TLine) > 0
                for k = 1:height(TLine)
                    i1 = find(busIds == round(TLine.fbus(k)), 1);
                    i2 = find(busIds == round(TLine.tbus(k)), 1);
                    if ~isempty(i1) && ~isempty(i2) && i1 ~= i2
                        sEdges(end+1) = i1; %#ok<AGROW>
                        tEdges(end+1) = i2; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(TIC) && height(TIC) > 0
                for k = 1:height(TIC)
                    i1 = find(busIds == round(TIC.fbus(k)), 1);
                    i2 = find(busIds == round(TIC.tbus(k)), 1);
                    if ~isempty(i1) && ~isempty(i2) && i1 ~= i2
                        sEdges(end+1) = i1; %#ok<AGROW>
                        tEdges(end+1) = i2; %#ok<AGROW>
                    end
                end
            end

            if strcmp(mode, 'circle') || isempty(sEdges)
                ang = linspace(0, 2*pi, nB + 1);
                ang(end) = [];
                coords = [cos(ang(:)), sin(ang(:))];
                return;
            end

            G = graph(sEdges, tEdges, [], nB);
            figTmp = figure('Visible', 'off');
            cleanupFig = onCleanup(@() close(figTmp));
            if isempty(cleanupFig), coords = zeros(nB, 2); return; end
            axTmp = axes(figTmp);
            if strcmp(mode, 'layered')
                h = plot(axTmp, G, 'Layout', 'layered');
            else
                h = plot(axTmp, G, 'Layout', 'force');
            end
            x = h.XData(:);
            y = h.YData(:);
            rx = max(x) - min(x); if rx == 0, rx = 1; end
            ry = max(y) - min(y); if ry == 0, ry = 1; end
            coords = [(x - min(x)) / rx, (y - min(y)) / ry];
        end

        function onCanvasClick(obj, evt)
            if isempty(obj.nodeCoords) || size(obj.nodeCoords, 1) == 0
                return;
            end
            pt = evt.IntersectionPoint(1:2);
            d2 = sum((obj.nodeCoords - pt).^2, 2);
            [~, idx] = min(d2);
            if idx >= 1 && idx <= numel(obj.ddBusInspect.Items)
                obj.ddBusInspect.Value = obj.ddBusInspect.Items{idx};
                obj.updateInspector();
            end
        end

        function updateInspector(obj)
            sess = obj.app.session;
            TBus = sess.tables.Bus;
            if isempty(TBus) || height(TBus) == 0 || strcmp(obj.ddBusInspect.Value, '-')
                return;
            end
            tok = regexp(char(obj.ddBusInspect.Value), '(\d+)', 'tokens', 'once');
            if isempty(tok), return; end
            bId = str2double(tok{1});
            rIdx = find(round(TBus.bus_i) == bId, 1);
            if isempty(rIdx), return; end

            isDC = startsWith(strtrim(TBus.Complexity_bus(rIdx)), "DC", 'IgnoreCase', true);
            domStr = 'AC';
            if isDC, domStr = 'DC'; end

            lines = {};
            lines{end+1} = sprintf('=== BUS %d (%s) ===', bId, domStr);
            lines{end+1} = sprintf('MATPOWER Type : %d', round(TBus.type(rIdx)));
            lines{end+1} = sprintf('Voltage Vm    : %.5f pu', TBus.Vm(rIdx));
            lines{end+1} = sprintf('Angle Va      : %.4f deg', TBus.Va(rIdx));
            lines{end+1} = sprintf('Load Pd / Qd  : %.3f MW / %.3f MVAr', TBus.Pd(rIdx), TBus.Qd(rIdx));
            lines{end+1} = sprintf('Bus Model     : %s', char(TBus.Complexity_bus(rIdx)));
            lines{end+1} = sprintf('Load Model    : %s', char(TBus.Complexity_load(rIdx)));
            lines{end+1} = '';

            % Connected Generators
            TGen = sess.tables.Gen;
            if ~isempty(TGen) && height(TGen) > 0
                gMask = find(round(TGen.bus) == bId);
                if ~isempty(gMask)
                    lines{end+1} = '--- Connected Units ---';
                    for g = gMask(:)'
                        lines{end+1} = sprintf('  * %s (mBase=%.1f MVA): Pg=%.2f MW, Qg=%.2f MVAr', ...
                            char(TGen.type(g)), TGen.mBase(g), TGen.Pg(g), TGen.Qg(g)); %#ok<AGROW>
                    end
                    lines{end+1} = '';
                end
            end

            % Connected Interlinking Converters
            TIC = sess.tables.IntConv;
            if ~isempty(TIC) && height(TIC) > 0
                icMask = find(round(TIC.fbus) == bId | round(TIC.tbus) == bId);
                if ~isempty(icMask)
                    lines{end+1} = '--- Interlinking Converters ---';
                    for c = icMask(:)'
                        lines{end+1} = sprintf('  * %s (%d -> %d, status=%d)', ...
                            char(TIC.type(c)), round(TIC.fbus(c)), round(TIC.tbus(c)), round(TIC.status(c))); %#ok<AGROW>
                    end
                    lines{end+1} = '';
                end
            end

            % Incident Branch Flows
            if ~isempty(sess.pfResult) && isfield(sess.pfResult, 'branchFlows')
                bf = sess.pfResult.branchFlows;
                bMask = find(round(bf.fbus) == bId | round(bf.tbus) == bId);
                if ~isempty(bMask)
                    lines{end+1} = '--- Solved Branch Flows ---';
                    for b = bMask(:)'
                        lines{end+1} = sprintf('  %d -> %d : Pf=%.2f MW, Qf=%.2f MVAr (Ploss=%.3f MW)', ...
                            round(bf.fbus(b)), round(bf.tbus(b)), bf.Pf_MW(b), bf.Qf_MVAr(b), bf.Ploss_MW(b)); %#ok<AGROW>
                    end
                end
            end

            obj.txtInspect.Value = lines;
        end

        function dlg = openDeviceForSelectedBus(obj)
            sess = obj.app.session;
            uName = '';
            tok = regexp(char(obj.ddBusInspect.Value), '(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                bId = str2double(tok{1});
                if ~isempty(sess.tables.Gen) && height(sess.tables.Gen) > 0
                    idx = find(round(sess.tables.Gen.bus) == bId, 1);
                    if ~isempty(idx)
                        uName = char(sess.tables.Gen.type(idx));
                    end
                end
                if isempty(uName) && ~isempty(sess.tables.IntConv) && height(sess.tables.IntConv) > 0
                    idx = find(round(sess.tables.IntConv.fbus) == bId | round(sess.tables.IntConv.tbus) == bId, 1);
                    if ~isempty(idx)
                        uName = char(sess.tables.IntConv.type(idx));
                    end
                end
            end
            dlg = DeviceDialog(obj.app, uName);
        end
    end
end
