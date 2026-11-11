classdef VariablesDialog < handle
    %VARIABLESDIALOG Symbolic $variable & sweep-range manager for VFlexPApp.

    properties
        app
        fig
        tbl
    end

    methods
        function obj = VariablesDialog(app)
            obj.app = app;
            obj.buildUI();
            obj.refresh();
        end

        function buildUI(obj)
            obj.fig = uifigure( ...
                'Name', 'VFlexP — Symbolic Variables ($name)', ...
                'Position', [240, 200, 660, 400], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [3, 1]);
            gl.RowHeight = {26, '1x', 38};
            gl.Padding = [12, 12, 12, 12];
            gl.RowSpacing = 8;

            uilabel(gl, ...
                'Text', 'Define symbolic $variables referenced in workbook cells (e.g. $P_load) or .dyr files (e.g. $H_vsm):', ...
                'FontColor', [0.25, 0.30, 0.38]);

            obj.tbl = uitable(gl, ...
                'ColumnName', {'Parameter', 'Value', 'Sweep From', 'Sweep To', 'Steps'}, ...
                'ColumnEditable', [true, true, true, true, true], ...
                'ColumnFormat', {'char', 'numeric', 'numeric', 'numeric', 'numeric'}, ...
                'RowName', 'numbered');

            bar = uigridlayout(gl, [1, 5]);
            bar.ColumnWidth = {115, 115, '1x', 110, 100};
            bar.Padding = [0, 0, 0, 0];

            uibutton(bar, 'Text', '+ Add Variable', ...
                'ButtonPushedFcn', @(~,~) obj.onAddRow());
            uibutton(bar, 'Text', '- Remove Last', ...
                'ButtonPushedFcn', @(~,~) obj.onRemoveRow());
            uilabel(bar, 'Text', '');
            uibutton(bar, 'Text', 'Apply', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.apply());
            uibutton(bar, 'Text', 'Close', ...
                'ButtonPushedFcn', @(~,~) obj.close());
        end

        function refresh(obj)
            if isempty(obj.fig) || ~isvalid(obj.fig)
                return;
            end
            sess = obj.app.session;
            vr = sess.variableRanges;
            if isempty(vr) || height(vr) == 0
                fn = fieldnames(sess.variables);
                if isempty(fn)
                    data = cell(0, 5);
                else
                    data = cell(numel(fn), 5);
                    for k = 1:numel(fn)
                        v = double(sess.variables.(fn{k}));
                        data(k, :) = {fn{k}, v, 0.8*v, 1.2*v, 5};
                    end
                end
            else
                n = height(vr);
                data = cell(n, 5);
                for k = 1:n
                    data(k, :) = {char(vr.Parameter(k)), double(vr.Value(k)), ...
                        double(vr.From(k)), double(vr.To(k)), round(double(vr.N_steps(k)))};
                end
            end
            obj.tbl.Data = data;
        end

        function onAddRow(obj)
            d = obj.tbl.Data;
            idx = size(d, 1) + 1;
            d(end+1, :) = {sprintf('var_%d', idx), 1.0, 0.8, 1.2, 5};
            obj.tbl.Data = d;
        end

        function onRemoveRow(obj)
            d = obj.tbl.Data;
            if ~isempty(d)
                d(end, :) = [];
                obj.tbl.Data = d;
            end
        end

        function apply(obj)
            sess = obj.app.session;
            d = obj.tbl.Data;
            n = size(d, 1);
            newVars = struct();
            pCol = strings(0, 1);
            vCol = zeros(0, 1);
            fCol = zeros(0, 1);
            tCol = zeros(0, 1);
            sCol = zeros(0, 1);
            for k = 1:n
                nm = regexprep(strtrim(char(string(d{k, 1}))), '^\$', '');
                if isempty(nm) || ~isvarname(nm)
                    continue;
                end
                vl = double(d{k, 2});
                fr = double(d{k, 3});
                to = double(d{k, 4});
                st = max(2, round(double(d{k, 5})));
                newVars.(nm) = vl;
                pCol(end+1, 1) = string(nm); %#ok<AGROW>
                vCol(end+1, 1) = vl; %#ok<AGROW>
                fCol(end+1, 1) = fr; %#ok<AGROW>
                tCol(end+1, 1) = to; %#ok<AGROW>
                sCol(end+1, 1) = st; %#ok<AGROW>
            end
            sess.variables = newVars;
            sess.variableRanges = table(pCol, vCol, fCol, tCol, sCol, ...
                'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
            for k = 1:numel(pCol)
                sess.setVariable(pCol(k), vCol(k));
            end
            sess.log('Updated %d symbolic $variables.', numel(pCol));
            obj.app.notifyVariablesChanged();
        end

        function close(obj)
            if ~isempty(obj.fig) && isvalid(obj.fig)
                delete(obj.fig);
            end
        end
    end
end
