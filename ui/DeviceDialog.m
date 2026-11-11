classdef DeviceDialog < handle
    %DEVICEDIALOG Per-unit submodel complexity & parameter override editor.

    properties
        app
        fig
        manifest
        ddUnit
        lblInfo
        ddMode1
        lblMode1
        ddMode2
        lblMode2
        tblParams
        mode1Field char = ''
        mode2Field char = ''
    end

    methods
        function obj = DeviceDialog(app, initialUnit)
            obj.app = app;
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            mfPath = fullfile(repoRoot, 'Simulink', 'lib', 'vflexp_lib_manifest.mat');
            S = load(mfPath, 'manifest');
            obj.manifest = S.manifest;

            if nargin < 2 || isempty(initialUnit)
                initialUnit = obj.pickDefaultUnit();
            end
            obj.buildUI(initialUnit);
            obj.refreshUnitView();
        end

        function uName = pickDefaultUnit(obj)
            sess = obj.app.session;
            allUnits = fieldnames(obj.manifest.units);
            if ~isempty(sess.tables.Gen) && height(sess.tables.Gen) > 0
                cand = char(sess.tables.Gen.type(1));
                if ismember(cand, allUnits)
                    uName = cand;
                    return;
                end
            end
            uName = allUnits{1};
        end

        function buildUI(obj, initialUnit)
            obj.fig = uifigure( ...
                'Name', 'VFlexP — Device Parameters & Complexity', ...
                'Position', [240, 160, 680, 520], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [4, 1]);
            gl.RowHeight = {36, 42, '1x', 38};
            gl.Padding = [12, 12, 12, 12];
            gl.RowSpacing = 8;

            % Row 1: Unit selector
            topBar = uigridlayout(gl, [1, 3]);
            topBar.ColumnWidth = {95, 220, '1x'};
            topBar.Padding = [0, 0, 0, 0];

            uilabel(topBar, 'Text', 'Unit model:', 'FontWeight', 'bold');
            allUnits = fieldnames(obj.manifest.units);
            if ~ismember(initialUnit, allUnits), initialUnit = allUnits{1}; end
            obj.ddUnit = uidropdown(topBar, ...
                'Items', allUnits, ...
                'Value', initialUnit, ...
                'ValueChangedFcn', @(~,~) obj.refreshUnitView());
            obj.lblInfo = uilabel(topBar, 'Text', '', 'FontColor', [0.25, 0.35, 0.48]);

            % Row 2: Submodel complexity selectors (shown dynamically per unit)
            modeBar = uigridlayout(gl, [1, 4]);
            modeBar.ColumnWidth = {140, 190, 140, '1x'};
            modeBar.Padding = [0, 0, 0, 0];
            obj.lblMode1 = uilabel(modeBar, 'Text', 'Submodel 1:');
            obj.ddMode1  = uidropdown(modeBar, 'Items', {'-'}, 'Value', '-');
            obj.lblMode2 = uilabel(modeBar, 'Text', 'Submodel 2:');
            obj.ddMode2  = uidropdown(modeBar, 'Items', {'-'}, 'Value', '-');

            % Row 3: Numeric parameters table
            obj.tblParams = uitable(gl, ...
                'ColumnName', {'Parameter', 'Default', 'Current Override / Value'}, ...
                'ColumnEditable', [false, false, true], ...
                'ColumnFormat', {'char', 'char', 'char'}, ...
                'RowName', 'numbered');

            % Row 4: Action buttons
            botBar = uigridlayout(gl, [1, 4]);
            botBar.ColumnWidth = {140, '1x', 110, 100};
            botBar.Padding = [0, 0, 0, 0];

            uibutton(botBar, 'Text', 'Reset Unit Overrides', ...
                'ButtonPushedFcn', @(~,~) obj.onResetUnit());
            uilabel(botBar, 'Text', '');
            uibutton(botBar, 'Text', 'Apply', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.apply());
            uibutton(botBar, 'Text', 'Close', ...
                'ButtonPushedFcn', @(~,~) obj.close());
        end

        function refreshUnitView(obj)
            uName = char(obj.ddUnit.Value);
            sess  = obj.app.session;
            defP  = obj.manifest.units.(uName).defaultPrm;
            ovP   = struct();
            if isfield(sess.unitOverrides, uName) && isstruct(sess.unitOverrides.(uName))
                ovP = sess.unitOverrides.(uName);
            end

            % Count active instances in current case
            nInst = 0;
            if ~isempty(sess.tables.Gen) && height(sess.tables.Gen) > 0
                nInst = nInst + sum(sess.tables.Gen.type == string(uName) & sess.tables.Gen.status > 0);
            end
            if ~isempty(sess.tables.IntConv) && height(sess.tables.IntConv) > 0
                nInst = nInst + sum(sess.tables.IntConv.type == string(uName) & sess.tables.IntConv.status > 0);
            end
            obj.lblInfo.Text = sprintf('Active instances in loaded case: %d', nInst);

            % Configure submodel complexity dropdowns
            obj.mode1Field = '';
            obj.mode2Field = '';
            obj.lblMode1.Visible = 'off'; obj.ddMode1.Visible = 'off';
            obj.lblMode2.Visible = 'off'; obj.ddMode2.Visible = 'off';

            if ismember(uName, {'MGES', 'SG'})
                obj.mode1Field = 'dynamic_stator';
                opts = {'order 8', 'order 4', 'order 2'};
                cur  = obj.resolveModeValue(uName, 'dynamic_stator', defP, ovP, 'order 8');
                obj.lblMode1.Text = 'Stator model:';
                obj.ddMode1.Items = opts;
                if ismember(cur, opts), obj.ddMode1.Value = cur; end
                obj.lblMode1.Visible = 'on'; obj.ddMode1.Visible = 'on';
            elseif ismember(uName, {'GFM', 'VSM', 'intconv_GFM', 'intconv_VSM'})
                obj.mode1Field = 'dynamic_VCC';
                if strcmp(uName, 'intconv_VSM')
                    opts = {'Static', 'Dynamic', 'Dynamic with current limitation', 'Dynamic with P limitation'};
                else
                    opts = {'Static', 'Dynamic', 'Dynamic with current limitation'};
                end
                cur = obj.resolveModeValue(uName, 'dynamic_VCC', defP, ovP, 'Static');
                obj.lblMode1.Text = 'Voltage/Current Ctrl:';
                obj.ddMode1.Items = opts;
                if ismember(cur, opts), obj.ddMode1.Value = cur; end
                obj.lblMode1.Visible = 'on'; obj.ddMode1.Visible = 'on';
            elseif ismember(uName, {'GFL', 'intconv_GFL'})
                obj.mode1Field = 'dynamic_CC';
                obj.mode2Field = 'dynamic_PLL';
                opts = {'Dynamic', 'Static'};
                cur1 = obj.resolveModeValue(uName, 'dynamic_CC', defP, ovP, 'Dynamic');
                cur2 = obj.resolveModeValue(uName, 'dynamic_PLL', defP, ovP, 'Dynamic');
                obj.lblMode1.Text = 'Current Ctrl (CC):';
                obj.ddMode1.Items = opts;
                if ismember(cur1, opts), obj.ddMode1.Value = cur1; end
                obj.lblMode2.Text = 'PLL model:';
                obj.ddMode2.Items = opts;
                if ismember(cur2, opts), obj.ddMode2.Value = cur2; end
                obj.lblMode1.Visible = 'on'; obj.ddMode1.Visible = 'on';
                obj.lblMode2.Visible = 'on'; obj.ddMode2.Visible = 'on';
            end

            % Populate numeric parameter rows
            fNames = fieldnames(defP);
            rows = cell(0, 3);
            for k = 1:numel(fNames)
                fn = fNames{k};
                if ismember(fn, {'dynamic_stator','dynamic_VCC','dynamic_CC','dynamic_PLL', ...
                                 'SGs_nodes','SCs_nodes','DFIGs_nodes'})
                    continue;
                end
                dv = defP.(fn);
                if isfield(ovP, fn)
                    cv = ovP.(fn);
                else
                    cv = dv;
                end
                rows(end+1, :) = {fn, mat2str(dv, 6), mat2str(cv, 6)}; %#ok<AGROW>
            end
            obj.tblParams.Data = rows;
        end

        function cur = resolveModeValue(~, ~, fld, defP, ovP, fallback)
            if isfield(ovP, fld)
                cur = char(string(ovP.(fld)));
            elseif isfield(defP, fld)
                cur = char(string(defP.(fld)));
            else
                cur = fallback;
            end
        end

        function onResetUnit(obj)
            uName = char(obj.ddUnit.Value);
            if isfield(obj.app.session.unitOverrides, uName)
                obj.app.session.unitOverrides = rmfield(obj.app.session.unitOverrides, uName);
            end
            obj.app.session.invalidateDynamicResults();
            obj.app.session.log('Cleared overrides for unit "%s".', uName);
            obj.refreshUnitView();
            obj.app.notifyComplexityChanged();
        end

        function apply(obj)
            uName = char(obj.ddUnit.Value);
            sess  = obj.app.session;
            defP  = obj.manifest.units.(uName).defaultPrm;
            ov    = struct();
            if isfield(sess.unitOverrides, uName) && isstruct(sess.unitOverrides.(uName))
                ov = sess.unitOverrides.(uName);
            end

            if ~isempty(obj.mode1Field)
                ov.(obj.mode1Field) = char(obj.ddMode1.Value);
            end
            if ~isempty(obj.mode2Field)
                ov.(obj.mode2Field) = char(obj.ddMode2.Value);
            end

            d = obj.tblParams.Data;
            for k = 1:size(d, 1)
                fn = d{k, 1};
                defStr = strtrim(d{k, 2});
                curStr = strtrim(d{k, 3});
                if ~strcmp(curStr, defStr)
                    val = str2num(curStr); %#ok<ST2NM>
                    if ~isempty(val)
                        ov.(fn) = double(val);
                    end
                elseif isfield(ov, fn) && isfield(defP, fn)
                    ov = rmfield(ov, fn);
                end
            end

            sess.unitOverrides.(uName) = ov;
            sess.invalidateDynamicResults();
            sess.log('Applied complexity/parameter overrides to unit "%s".', uName);
            obj.app.notifyComplexityChanged();
        end

        function close(obj)
            if ~isempty(obj.fig) && isvalid(obj.fig)
                delete(obj.fig);
            end
        end
    end
end
