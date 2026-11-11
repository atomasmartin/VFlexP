classdef StaticTab < handle
    %STATICTAB Editable network sheets, Power Flow & OPF tab for VFlexPApp.

    properties
        app
        parent
        edtSb
        lblSummary
        subTabs
        tblBus
        tblGen
        tblLine
        tblIntConv
        tblFb
        tblFlows
        selectedRow struct = struct('Bus', 0, 'Gen', 0, 'Line', 0, 'IntConv', 0, 'Fb', 0)
    end

    methods
        function obj = StaticTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [3, 1]);
            gl.RowHeight = {40, 28, '1x'};
            gl.Padding = [10, 10, 10, 10];
            gl.RowSpacing = 6;

            % Action toolbar
            bar = uigridlayout(gl, [1, 8]);
            bar.ColumnWidth = {115, 115, 72, 85, 105, 110, 105, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uibutton(bar, 'Text', 'Run PF', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.onRunPF());
            uibutton(bar, 'Text', 'Run OPF', ...
                'BackgroundColor', [0.16, 0.56, 0.36], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.onRunOPF());

            uilabel(bar, 'Text', 'Sb (MVA):', 'HorizontalAlignment', 'right');
            obj.edtSb = uieditfield(bar, 'numeric', ...
                'Value', 100, ...
                'Limits', [1e-3, 1e6], ...
                'ValueChangedFcn', @(~,~) obj.onSbChanged());

            uibutton(bar, 'Text', '+ Add Row', ...
                'ButtonPushedFcn', @(~,~) obj.onAddRow());
            uibutton(bar, 'Text', '- Delete Row', ...
                'ButtonPushedFcn', @(~,~) obj.onDeleteRow());
            uibutton(bar, 'Text', 'Settings...', ...
                'ButtonPushedFcn', @(~,~) StaticSettingsDialog(obj.app));
            uilabel(bar, 'Text', '');

            % KPI Summary Banner
            obj.lblSummary = uilabel(gl, ...
                'Text', 'Load a case and click "Run PF" or "Run OPF" to solve the steady-state operating point.', ...
                'FontWeight', 'bold', ...
                'FontColor', [0.18, 0.25, 0.35]);

            % Sub-tabs for workbook tables + solved branch flows
            obj.subTabs = uitabgroup(gl);

            tabBus = uitab(obj.subTabs, 'Title', 'Buses', 'Tag', 'Bus');
            g1 = uigridlayout(tabBus, [1, 1]); g1.Padding = [4, 4, 4, 4];
            obj.tblBus = uitable(g1, ...
                'CellEditCallback', @(~, evt) obj.onTableEdited('Bus', evt), ...
                'CellSelectionCallback', @(~, evt) obj.onRowSelected('Bus', evt));

            tabGen = uitab(obj.subTabs, 'Title', 'Generators / Units', 'Tag', 'Gen');
            g2 = uigridlayout(tabGen, [1, 1]); g2.Padding = [4, 4, 4, 4];
            obj.tblGen = uitable(g2, ...
                'CellEditCallback', @(~, evt) obj.onTableEdited('Gen', evt), ...
                'CellSelectionCallback', @(~, evt) obj.onRowSelected('Gen', evt));

            tabLine = uitab(obj.subTabs, 'Title', 'Lines & Transformers', 'Tag', 'Line');
            g3 = uigridlayout(tabLine, [1, 1]); g3.Padding = [4, 4, 4, 4];
            obj.tblLine = uitable(g3, ...
                'CellEditCallback', @(~, evt) obj.onTableEdited('Line', evt), ...
                'CellSelectionCallback', @(~, evt) obj.onRowSelected('Line', evt));

            tabIC = uitab(obj.subTabs, 'Title', 'Interlinking Converters', 'Tag', 'IntConv');
            g4 = uigridlayout(tabIC, [1, 1]); g4.Padding = [4, 4, 4, 4];
            obj.tblIntConv = uitable(g4, ...
                'CellEditCallback', @(~, evt) obj.onTableEdited('IntConv', evt), ...
                'CellSelectionCallback', @(~, evt) obj.onRowSelected('IntConv', evt));

            tabFb = uitab(obj.subTabs, 'Title', 'Areas (Fb)', 'Tag', 'Fb');
            g5 = uigridlayout(tabFb, [1, 1]); g5.Padding = [4, 4, 4, 4];
            obj.tblFb = uitable(g5, ...
                'CellEditCallback', @(~, evt) obj.onTableEdited('Fb', evt), ...
                'CellSelectionCallback', @(~, evt) obj.onRowSelected('Fb', evt));

            tabFlw = uitab(obj.subTabs, 'Title', 'Solved Branch Flows', 'Tag', 'Flows');
            g6 = uigridlayout(tabFlw, [1, 1]); g6.Padding = [4, 4, 4, 4];
            obj.tblFlows = uitable(g6, 'ColumnEditable', false);
        end

        function refresh(obj)
            sess = obj.app.session;
            obj.edtSb.Value = double(sess.tables.Sb);

            obj.populateEditableTable(obj.tblBus, sess.rawTables.Bus, VFlexPSession.busColumnNames());
            obj.populateEditableTable(obj.tblGen, sess.rawTables.Gen, VFlexPSession.genColumnNames());
            obj.populateEditableTable(obj.tblLine, sess.rawTables.Line, VFlexPSession.lineColumnNames());
            obj.populateEditableTable(obj.tblIntConv, sess.rawTables.IntConv, VFlexPSession.intConvColumnNames());
            obj.populateEditableTable(obj.tblFb, sess.rawTables.Fb, {'Area', 'Fb'});

            obj.refreshFlows();
        end

        function refreshFlows(obj)
            pf = obj.app.session.pfResult;
            if isempty(pf)
                obj.lblSummary.Text = 'Power Flow not yet solved for current tables.';
                obj.lblSummary.FontColor = [0.35, 0.38, 0.44];
                obj.tblFlows.Data = cell(0, 10);
                obj.tblFlows.ColumnName = {'fbus','tbus','Pf (MW)','Qf (MVAr)','Pt (MW)','Qt (MVAr)','Ploss (MW)','Loading (%)','DC','IntConv'};
                return;
            end

            s = pf.summary;
            if pf.converged
                statusTxt = sprintf('%s CONVERGED (%.3f s)', pf.kind, pf.elapsed);
                obj.lblSummary.FontColor = [0.08, 0.46, 0.24];
            else
                statusTxt = sprintf('%s DID NOT CONVERGE', pf.kind);
                obj.lblSummary.FontColor = [0.75, 0.15, 0.15];
            end
            costTxt = '';
            if isfield(s, 'cost')
                costTxt = sprintf('  |  Cost = %.2f $/h', s.cost);
            end
            obj.lblSummary.Text = sprintf( ...
                '%s  |  Pgen = %.2f MW, Qgen = %.2f MVAr  |  Pload = %.2f MW, Qload = %.2f MVAr  |  Ploss = %.3f MW  |  Vm in [%.4f, %.4f] pu%s', ...
                statusTxt, s.Pgen_MW, s.Qgen_MVAr, s.Pload_MW, s.Qload_MVAr, s.Ploss_MW, s.Vm_min, s.Vm_max, costTxt);

            bf = pf.branchFlows;
            n = height(bf);
            rows = cell(n, 10);
            for r = 1:n
                rows(r, :) = {bf.fbus(r), bf.tbus(r), bf.Pf_MW(r), bf.Qf_MVAr(r), ...
                    bf.Pt_MW(r), bf.Qt_MVAr(r), bf.Ploss_MW(r), bf.loading_pct(r), ...
                    logical(bf.isDC(r)), logical(bf.isIntConv(r))};
            end
            obj.tblFlows.ColumnName = {'fbus','tbus','Pf (MW)','Qf (MVAr)','Pt (MW)','Qt (MVAr)','Ploss (MW)','Loading (%)','DC','IntConv'};
            obj.tblFlows.Data = rows;
        end

        function populateEditableTable(~, uiTbl, rawTbl, colNames)
            nCols = numel(colNames);
            if isempty(rawTbl) || height(rawTbl) == 0
                uiTbl.ColumnName = colNames;
                uiTbl.ColumnEditable = true(1, nCols);
                uiTbl.Data = cell(0, nCols);
                return;
            end
            nRows = height(rawTbl);
            data = cell(nRows, nCols);
            for c = 1:nCols
                cn = colNames{c};
                if ismember(cn, rawTbl.Properties.VariableNames)
                    col = rawTbl.(cn);
                    for r = 1:nRows
                        if iscell(col)
                            v = col{r};
                        else
                            v = col(r);
                        end
                        if isnumeric(v) || islogical(v)
                            data{r, c} = num2str(double(v), '%.12g');
                        else
                            sv = string(v);
                            if ismissing(sv)
                                data{r, c} = '';
                            else
                                data{r, c} = char(sv);
                            end
                        end
                    end
                else
                    data(:, c) = {''};
                end
            end
            uiTbl.ColumnName = colNames;
            uiTbl.ColumnEditable = true(1, nCols);
            uiTbl.Data = data;
        end

        function onRowSelected(obj, sheetTag, evt)
            if ~isempty(evt.Indices)
                obj.selectedRow.(sheetTag) = evt.Indices(1, 1);
            end
        end

        function onSbChanged(obj)
            sess = obj.app.session;
            sess.rawTables.Sb = double(obj.edtSb.Value);
            sess.evalTables();
            sess.invalidateDynamicResults();
            sess.log('Updated system base Sb = %.4g MVA.', sess.tables.Sb);
            obj.app.refreshLogStrip();
        end

        function onTableEdited(obj, sheetName, evt)
            sess = obj.app.session;
            r = evt.Indices(1);
            c = evt.Indices(2);
            uiTbl = obj.getTableWidget(sheetName);
            colName = uiTbl.ColumnName{c};
            newStr = strtrim(char(string(evt.NewData)));

            RT = sess.rawTables.(sheetName);
            if ~ismember(colName, RT.Properties.VariableNames)
                RT.(colName) = strings(height(RT), 1);
            end
            col = RT.(colName);
            numVal = str2double(newStr);
            if isnumeric(col) && ~isnan(numVal) && ~startsWith(newStr, '$')
                RT.(colName)(r) = numVal;
            else
                sCol = string(col);
                sCol(r) = string(newStr);
                RT.(colName) = sCol;
            end
            sess.rawTables.(sheetName) = RT;
            sess.evalTables();
            sess.invalidateDynamicResults();
            sess.log('Edited %s(%d, %s) = "%s".', sheetName, r, colName, newStr);
            obj.app.notifyComplexityChanged();
        end

        function onAddRow(obj)
            tag = char(obj.subTabs.SelectedTab.Tag);
            if strcmp(tag, 'Flows')
                return;
            end
            sess = obj.app.session;
            switch tag
                case 'Bus'
                    T = sess.tables.Bus;
                    nextId = 1;
                    if ~isempty(T) && height(T) > 0
                        nextId = max(T.bus_i) + 1;
                    end
                    row = {num2str(nextId), '1', '0', '0', '0', '0', '1', '1', '0', '230', '1', '1.1', '0.9', ...
                           'static omega', 'series RL static omega', '0', '0'};
                    cols = VFlexPSession.busColumnNames();
                case 'Gen'
                    bId = 1;
                    if ~isempty(sess.tables.Bus) && height(sess.tables.Bus) > 0
                        bId = sess.tables.Bus.bus_i(1);
                    end
                    cols = VFlexPSession.genColumnNames();
                    row = repmat({'0'}, 1, numel(cols));
                    row{1} = num2str(bId);
                    row{2} = 'VSM';
                    row{5} = '100';
                    row{6} = '-100';
                    row{7} = '1';
                    row{8} = num2str(sess.tables.Sb);
                    row{9} = '1';
                    row{10} = '100';
                case 'Line'
                    cols = VFlexPSession.lineColumnNames();
                    row = {'1', '2', '0.01', '0.1', '0', '250', '250', '250', '0', '0', '1', '-360', '360', 'static omega'};
                case 'IntConv'
                    cols = VFlexPSession.intConvColumnNames();
                    row = repmat({'0'}, 1, numel(cols));
                    row{1} = 'intconv_GFL';
                    row{2} = '1';
                    row{3} = '2';
                    row{4} = '0.005';
                    row{5} = '0.05';
                    row{12} = '1';
                    row{25} = '1.1';
                    row{26} = '0.9';
                    row{27} = '1';
                    row{29} = '1';
                    row{32} = '-360';
                    row{33} = '360';
                case 'Fb'
                    cols = {'Area', 'Fb'};
                    nextA = height(sess.tables.Fb) + 1;
                    row = {num2str(nextA), '50'};
                otherwise
                    return;
            end
            uiTbl = obj.getTableWidget(tag);
            d = uiTbl.Data;
            d(end+1, :) = row;
            uiTbl.Data = d;
            sess.rawTables.(tag) = cell2table(d, 'VariableNames', cols);
            sess.evalTables();
            sess.invalidateDynamicResults();
            sess.log('Added row %d to %s table.', size(d, 1), tag);
            obj.app.notifyComplexityChanged();
        end

        function onDeleteRow(obj)
            tag = char(obj.subTabs.SelectedTab.Tag);
            if strcmp(tag, 'Flows')
                return;
            end
            uiTbl = obj.getTableWidget(tag);
            d = uiTbl.Data;
            if isempty(d)
                return;
            end
            r = obj.selectedRow.(tag);
            if r < 1 || r > size(d, 1)
                r = size(d, 1);
            end
            d(r, :) = [];
            uiTbl.Data = d;
            sess = obj.app.session;
            sess.rawTables.(tag) = cell2table(d, 'VariableNames', uiTbl.ColumnName);
            sess.evalTables();
            sess.invalidateDynamicResults();
            sess.log('Deleted row %d from %s table.', r, tag);
            obj.app.notifyComplexityChanged();
        end

        function res = onRunPF(obj)
            res = obj.app.runPF();
        end

        function res = onRunOPF(obj)
            res = obj.app.runOPF();
        end

        function w = getTableWidget(obj, tag)
            switch tag
                case 'Bus',     w = obj.tblBus;
                case 'Gen',     w = obj.tblGen;
                case 'Line',    w = obj.tblLine;
                case 'IntConv', w = obj.tblIntConv;
                case 'Fb',      w = obj.tblFb;
                otherwise,      w = obj.tblFlows;
            end
        end
    end
end
