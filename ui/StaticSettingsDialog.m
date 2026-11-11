classdef StaticSettingsDialog < handle
    %STATICSETTINGSDIALOG Power-flow & OPF configuration dialog for VFlexPApp.

    properties
        app
        fig
        edtTol
        spnMaxIter
        chkEnforceQ
        chkWriteBranch
    end

    methods
        function obj = StaticSettingsDialog(app)
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            cfg = obj.app.session.config.powerflow;
            obj.fig = uifigure( ...
                'Name', 'Static Analysis Settings', ...
                'Position', [300, 260, 420, 250], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [5, 2]);
            gl.RowHeight = {32, 32, 30, 30, '1x'};
            gl.ColumnWidth = {190, '1x'};
            gl.Padding = [14, 14, 14, 14];
            gl.RowSpacing = 8;

            uilabel(gl, 'Text', 'Convergence tolerance:');
            obj.edtTol = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.tolerance), ...
                'ValueDisplayFormat', '%.2e', ...
                'Limits', [1e-14, 1e-1]);

            uilabel(gl, 'Text', 'Maximum iterations:');
            obj.spnMaxIter = uispinner(gl, ...
                'Value', round(double(cfg.max_iterations)), ...
                'Limits', [1, 500], ...
                'Step', 5);

            obj.chkEnforceQ = uicheckbox(gl, ...
                'Text', 'Enforce generator Q limits in PF', ...
                'Value', logical(cfg.enforce_q_limits));
            obj.chkEnforceQ.Layout.Column = [1, 2];

            obj.chkWriteBranch = uicheckbox(gl, ...
                'Text', 'Write solved flows back to raw tables', ...
                'Value', logical(cfg.write_branch));
            obj.chkWriteBranch.Layout.Column = [1, 2];

            bar = uigridlayout(gl, [1, 3]);
            bar.Layout.Row = 5;
            bar.Layout.Column = [1, 2];
            bar.ColumnWidth = {'1x', 100, 90};
            bar.Padding = [0, 4, 0, 0];

            uilabel(bar, 'Text', '');
            uibutton(bar, 'Text', 'Apply', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.applyAndClose());
            uibutton(bar, 'Text', 'Cancel', ...
                'ButtonPushedFcn', @(~,~) obj.close());
        end

        function apply(obj)
            sess = obj.app.session;
            sess.config.powerflow.tolerance = double(obj.edtTol.Value);
            sess.config.powerflow.max_iterations = round(double(obj.spnMaxIter.Value));
            sess.config.powerflow.enforce_q_limits = logical(obj.chkEnforceQ.Value);
            sess.config.powerflow.write_branch = logical(obj.chkWriteBranch.Value);
            sess.log('Updated Power Flow settings (tol=%.2e, max_iter=%d).', ...
                sess.config.powerflow.tolerance, sess.config.powerflow.max_iterations);
            obj.app.refreshLogStrip();
        end

        function applyAndClose(obj)
            obj.apply();
            obj.close();
        end

        function close(obj)
            if ~isempty(obj.fig) && isvalid(obj.fig)
                delete(obj.fig);
            end
        end
    end
end
