classdef ComplexityTab < handle
    %COMPLEXITYTAB Element-by-element & unit-submodel complexity manager tab.

    properties
        app
        parent
        lblSummary
        tblBusComp
        tblLineComp
        tblUnitComp
        manifest
    end

    methods
        function obj = ComplexityTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            S = load(fullfile(repoRoot, 'Simulink', 'lib', 'vflexp_lib_manifest.mat'), 'manifest');
            obj.manifest = S.manifest;
            obj.buildUI();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [3, 1]);
            gl.RowHeight = {40, 28, '1x'};
            gl.Padding = [10, 10, 10, 10];
            gl.RowSpacing = 6;

            % Toolbar
            bar = uigridlayout(gl, [1, 6]);
            bar.ColumnWidth = {165, 165, 160, 155, '1x', 95};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uibutton(bar, 'Text', 'Preset: All Dynamic (EMT)', ...
                'ButtonPushedFcn', @(~,~) obj.applyPreset('EMT'));
            uibutton(bar, 'Text', 'Preset: All Static (RMS)', ...
                'ButtonPushedFcn', @(~,~) obj.applyPreset('RMS'));
            uibutton(bar, 'Text', 'Verify & Count States', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.verifyAndCountStates());
            uibutton(bar, 'Text', 'Device Parameters...', ...
                'ButtonPushedFcn', @(~,~) DeviceDialog(obj.app));
            uilabel(bar, 'Text', '');
            uibutton(bar, 'Text', 'Refresh', ...
                'ButtonPushedFcn', @(~,~) obj.refresh());

            obj.lblSummary = uilabel(gl, ...
                'Text', 'Edit bus, load, branch, or unit submodel complexity below.', ...
                'FontWeight', 'bold', ...
                'FontColor', [0.18, 0.25, 0.35]);

            % 3-panel layout: Buses & Loads | Lines | Active Units
            body = uigridlayout(gl, [1, 3]);
            body.ColumnWidth = {'1.2x', '1x', '1.1x'};
            body.Padding = [0, 0, 0, 0];
            body.ColumnSpacing = 8;

            pBus = uipanel(body, 'Title', 'Buses & Shunt Loads');
            gb = uigridlayout(pBus, [1, 1]); gb.Padding = [4, 4, 4, 4];
            busPop = {'dynamic', 'static omega', 'static', 'only R', 'DC dynamic', 'DC static'};
            loadPop = {'series RL dynamic', 'series RL static omega', 'series RL static', ...
                       'PQ', 'dynamic PQ', 'parallel RL dynamic', 'parallel RL static omega', ...
                       'parallel RL static', 'constant current', ...
                       'DC ConstR', 'DC ConstP', 'DC ConstI', 'DC DynamicP'};
            obj.tblBusComp = uitable(gb, ...
                'ColumnName', {'Bus', 'Domain', 'Bus Complexity', 'Load Complexity'}, ...
                'ColumnEditable', [false, false, true, true], ...
                'ColumnFormat', {'numeric', 'char', busPop, loadPop}, ...
                'RowName', [], ...
                'CellEditCallback', @(~, evt) obj.onBusCompEdited(evt));

            pLine = uipanel(body, 'Title', 'AC / DC Lines & Transformers');
            gln = uigridlayout(pLine, [1, 1]); gln.Padding = [4, 4, 4, 4];
            linePop = {'dynamic', 'static omega', 'static', 'DC dynamic', 'DC static'};
            obj.tblLineComp = uitable(gln, ...
                'ColumnName', {'From', 'To', 'Status', 'Branch Complexity'}, ...
                'ColumnEditable', [false, false, false, true], ...
                'ColumnFormat', {'numeric', 'numeric', 'numeric', linePop}, ...
                'RowName', [], ...
                'CellEditCallback', @(~, evt) obj.onLineCompEdited(evt));

            pUnit = uipanel(body, 'Title', 'Active Generators & Converters');
            gu = uigridlayout(pUnit, [1, 1]); gu.Padding = [4, 4, 4, 4];
            unitPop = {'order 8', 'order 4', 'order 2', 'Dynamic', 'Static', ...
                       'Dynamic with current limitation', 'Dynamic with P limitation', 'Fixed (Standard)'};
            obj.tblUnitComp = uitable(gu, ...
                'ColumnName', {'Unit Type', 'Count', 'Submodel', 'Mode'}, ...
                'ColumnEditable', [false, false, false, true], ...
                'ColumnFormat', {'char', 'numeric', 'char', unitPop}, ...
                'RowName', [], ...
                'CellEditCallback', @(~, evt) obj.onUnitCompEdited(evt));
        end

        function refresh(obj)
            sess = obj.app.session;
            TBus = sess.tables.Bus;
            nB = height(TBus);
            bRows = cell(nB, 4);
            for r = 1:nB
                cb = char(strtrim(TBus.Complexity_bus(r)));
                cl = char(strtrim(TBus.Complexity_load(r)));
                dom = 'AC';
                if startsWith(cb, 'DC', 'IgnoreCase', true), dom = 'DC'; end
                bRows(r, :) = {TBus.bus_i(r), dom, cb, cl};
            end
            obj.tblBusComp.Data = bRows;

            TLine = sess.tables.Line;
            nL = height(TLine);
            lRows = cell(nL, 4);
            for r = 1:nL
                lRows(r, :) = {TLine.fbus(r), TLine.tbus(r), TLine.status(r), char(strtrim(TLine.Complexity(r)))};
            end
            obj.tblLineComp.Data = lRows;

            % Active unit types
            uRows = cell(0, 4);
            activeTypes = strings(0, 1);
            counts = zeros(0, 1);
            if ~isempty(sess.tables.Gen) && height(sess.tables.Gen) > 0
                onG = sess.tables.Gen(sess.tables.Gen.status > 0, :);
                uG = unique(onG.type, 'stable');
                for k = 1:numel(uG)
                    activeTypes(end+1, 1) = uG(k); %#ok<AGROW>
                    counts(end+1, 1) = sum(onG.type == uG(k)); %#ok<AGROW>
                end
            end
            if ~isempty(sess.tables.IntConv) && height(sess.tables.IntConv) > 0
                onC = sess.tables.IntConv(sess.tables.IntConv.status > 0, :);
                uC = unique(onC.type, 'stable');
                for k = 1:numel(uC)
                    activeTypes(end+1, 1) = uC(k); %#ok<AGROW>
                    counts(end+1, 1) = sum(onC.type == uC(k)); %#ok<AGROW>
                end
            end

            for k = 1:numel(activeTypes)
                uName = char(activeTypes(k));
                cnt = counts(k);
                if ~isfield(obj.manifest.units, uName), continue; end
                defP = obj.manifest.units.(uName).defaultPrm;
                ovP = struct();
                if isfield(sess.unitOverrides, uName)
                    ovP = sess.unitOverrides.(uName);
                end
                if ismember(uName, {'MGES', 'SG'})
                    cur = obj.getMode(defP, ovP, 'dynamic_stator', 'order 8');
                    uRows(end+1, :) = {uName, cnt, 'dynamic_stator', cur}; %#ok<AGROW>
                elseif ismember(uName, {'GFM', 'VSM', 'intconv_GFM', 'intconv_VSM'})
                    cur = obj.getMode(defP, ovP, 'dynamic_VCC', 'Static');
                    uRows(end+1, :) = {uName, cnt, 'dynamic_VCC', cur}; %#ok<AGROW>
                elseif ismember(uName, {'GFL', 'intconv_GFL'})
                    cur1 = obj.getMode(defP, ovP, 'dynamic_CC', 'Dynamic');
                    cur2 = obj.getMode(defP, ovP, 'dynamic_PLL', 'Dynamic');
                    uRows(end+1, :) = {uName, cnt, 'dynamic_CC', cur1}; %#ok<AGROW>
                    uRows(end+1, :) = {uName, cnt, 'dynamic_PLL', cur2}; %#ok<AGROW>
                else
                    uRows(end+1, :) = {uName, cnt, '-', 'Fixed (Standard)'}; %#ok<AGROW>
                end
            end
            obj.tblUnitComp.Data = uRows;

            if ~isempty(sess.initResult)
                obj.lblSummary.Text = sprintf( ...
                    'Initialized Dynamic State Count: %d states  |  AC Grid Mode: %d  |  DC Grid Mode: %d', ...
                    sess.initResult.nStates, sess.initResult.acGridMode, sess.initResult.dcGridMode);
            else
                nDynBus = sum(contains(lower(TBus.Complexity_bus), "dynamic"));
                nDynLine = sum(contains(lower(TLine.Complexity), "dynamic"));
                obj.lblSummary.Text = sprintf( ...
                    'Dynamic buses: %d/%d  |  Dynamic branches: %d/%d  |  Click "Verify & Count States" to compile state count.', ...
                    nDynBus, nB, nDynLine, nL);
            end
        end

        function cur = getMode(~, defP, ovP, fld, fallback)
            if isfield(ovP, fld)
                cur = char(string(ovP.(fld)));
            elseif isfield(defP, fld)
                cur = char(string(defP.(fld)));
            else
                cur = fallback;
            end
        end

        function onBusCompEdited(obj, evt)
            sess = obj.app.session;
            r = evt.Indices(1);
            c = evt.Indices(2);
            val = string(evt.NewData);
            if c == 3
                sess.tables.Bus.Complexity_bus(r) = val;
                sess.rawTables.Bus.Complexity_bus = string(sess.tables.Bus.Complexity_bus);
            elseif c == 4
                sess.tables.Bus.Complexity_load(r) = val;
                sess.rawTables.Bus.Complexity_load = string(sess.tables.Bus.Complexity_load);
            end
            sess.invalidateDynamicResults();
            sess.log('Updated Bus %d complexity -> (%s, %s).', ...
                sess.tables.Bus.bus_i(r), sess.tables.Bus.Complexity_bus(r), sess.tables.Bus.Complexity_load(r));
            obj.app.notifyComplexityChanged();
        end

        function onLineCompEdited(obj, evt)
            sess = obj.app.session;
            r = evt.Indices(1);
            val = string(evt.NewData);
            sess.tables.Line.Complexity(r) = val;
            sess.rawTables.Line.Complexity = string(sess.tables.Line.Complexity);
            sess.invalidateDynamicResults();
            sess.log('Updated Line (%d-%d) complexity -> %s.', ...
                sess.tables.Line.fbus(r), sess.tables.Line.tbus(r), val);
            obj.app.notifyComplexityChanged();
        end

        function onUnitCompEdited(obj, evt)
            sess = obj.app.session;
            r = evt.Indices(1);
            d = obj.tblUnitComp.Data;
            uName = d{r, 1};
            subFld = d{r, 3};
            newVal = char(string(evt.NewData));
            if strcmp(subFld, '-')
                obj.refresh();
                return;
            end
            if ~isfield(sess.unitOverrides, uName) || ~isstruct(sess.unitOverrides.(uName))
                sess.unitOverrides.(uName) = struct();
            end
            sess.unitOverrides.(uName).(subFld) = newVal;
            sess.invalidateDynamicResults();
            sess.log('Updated %s.%s -> "%s".', uName, subFld, newVal);
            obj.app.notifyComplexityChanged();
        end

        function applyPreset(obj, presetName)
            sess = obj.app.session;
            isEMT = strcmpi(presetName, 'EMT');
            for r = 1:height(sess.tables.Bus)
                isDC = startsWith(strtrim(sess.tables.Bus.Complexity_bus(r)), "DC", 'IgnoreCase', true);
                if isDC
                    curLd = lower(strtrim(sess.tables.Bus.Complexity_load(r)));
                    if isEMT
                        sess.tables.Bus.Complexity_bus(r) = "DC dynamic";
                        if contains(curLd, "constp") || contains(curLd, "dynamicp")
                            if sess.tables.Bus.Bs(r) >= 0.05
                                sess.tables.Bus.Complexity_load(r) = "DC DynamicP";
                            end
                        elseif curLd == "" || ~startsWith(curLd, "dc")
                            sess.tables.Bus.Complexity_load(r) = "DC ConstR";
                        end
                    else
                        sess.tables.Bus.Complexity_bus(r) = "DC static";
                        if contains(curLd, "dynamicp")
                            sess.tables.Bus.Complexity_load(r) = "DC ConstP";
                        elseif curLd == "" || ~startsWith(curLd, "dc")
                            sess.tables.Bus.Complexity_load(r) = "DC ConstR";
                        end
                    end
                else
                    if isEMT
                        sess.tables.Bus.Complexity_bus(r) = "dynamic";
                        if sess.tables.Bus.Qd(r) > 0
                            sess.tables.Bus.Complexity_load(r) = "series RL dynamic";
                        else
                            sess.tables.Bus.Complexity_load(r) = "series RL static omega";
                        end
                    else
                        sess.tables.Bus.Complexity_bus(r) = "static omega";
                        sess.tables.Bus.Complexity_load(r) = "series RL static omega";
                    end
                end
            end
            sess.rawTables.Bus.Complexity_bus = string(sess.tables.Bus.Complexity_bus);
            sess.rawTables.Bus.Complexity_load = string(sess.tables.Bus.Complexity_load);

            for r = 1:height(sess.tables.Line)
                isDC = startsWith(strtrim(sess.tables.Line.Complexity(r)), "DC", 'IgnoreCase', true);
                if isDC
                    if isEMT && sess.tables.Line.x(r) > 0
                        sess.tables.Line.Complexity(r) = "DC dynamic";
                    else
                        sess.tables.Line.Complexity(r) = "DC static";
                    end
                else
                    if isEMT
                        sess.tables.Line.Complexity(r) = "dynamic";
                    else
                        sess.tables.Line.Complexity(r) = "static omega";
                    end
                end
            end
            if ~isempty(sess.tables.Line) && height(sess.tables.Line) > 0
                sess.rawTables.Line.Complexity = string(sess.tables.Line.Complexity);
            end

            % Unit submodels
            for u = {'MGES', 'SG'}
                uName = u{1};
                if ~isfield(sess.unitOverrides, uName), sess.unitOverrides.(uName) = struct(); end
                if isEMT
                    sess.unitOverrides.(uName).dynamic_stator = 'order 8';
                else
                    sess.unitOverrides.(uName).dynamic_stator = 'order 4';
                end
            end
            for u = {'GFM', 'VSM', 'intconv_GFM', 'intconv_VSM'}
                uName = u{1};
                if ~isfield(sess.unitOverrides, uName), sess.unitOverrides.(uName) = struct(); end
                if isEMT
                    sess.unitOverrides.(uName).dynamic_VCC = 'Dynamic';
                else
                    sess.unitOverrides.(uName).dynamic_VCC = 'Static';
                end
            end
            for u = {'GFL', 'intconv_GFL'}
                uName = u{1};
                if ~isfield(sess.unitOverrides, uName), sess.unitOverrides.(uName) = struct(); end
                sess.unitOverrides.(uName).dynamic_CC = 'Dynamic';
                sess.unitOverrides.(uName).dynamic_PLL = 'Dynamic';
            end

            sess.invalidateDynamicResults();
            sess.log('Applied "%s" complexity preset across network and unit submodels.', upper(presetName));
            obj.app.notifyComplexityChanged();
        end

        function initRes = verifyAndCountStates(obj)
            sess = obj.app.session;
            if isempty(sess.pfResult) || ~sess.pfResult.converged
                vflexp_run_pf(sess);
            end
            initRes = vflexp_build_and_init(sess);
            obj.refresh();
            obj.app.refreshLogStrip();
        end
    end
end
