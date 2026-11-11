classdef DynamicSettingsDialog < handle
    %DYNAMICSETTINGSDIALOG Time-domain solver & secondary/tertiary control dialog.

    properties
        app
        fig
        edtTend
        ddMethod
        edtRtol
        edtAtol
        edtMaxStep
        ddSC
        ddTC
    end

    methods
        function obj = DynamicSettingsDialog(app)
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            cfg = obj.app.session.config.sim;
            obj.fig = uifigure( ...
                'Name', 'Dynamic Simulation Settings', ...
                'Position', [290, 210, 450, 350], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [8, 2]);
            gl.RowHeight = {30, 30, 30, 30, 30, 30, 30, '1x'};
            gl.ColumnWidth = {200, '1x'};
            gl.Padding = [14, 14, 14, 14];
            gl.RowSpacing = 8;

            uilabel(gl, 'Text', 'Stop time t_end (s):');
            obj.edtTend = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.t_end), 'Limits', [1e-3, 1e5]);

            uilabel(gl, 'Text', 'ODE/DAE solver:');
            solvers = {'ode15s', 'ode23tb', 'ode23s', 'ode23t', 'ode45'};
            curSol = char(string(cfg.method));
            if ~ismember(curSol, solvers), solvers{end+1} = curSol; end
            obj.ddMethod = uidropdown(gl, 'Items', solvers, 'Value', curSol);

            uilabel(gl, 'Text', 'Relative tolerance (RelTol):');
            obj.edtRtol = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.rtol), 'ValueDisplayFormat', '%.2e', 'Limits', [1e-12, 1e-1]);

            uilabel(gl, 'Text', 'Absolute tolerance (AbsTol):');
            obj.edtAtol = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.atol), 'ValueDisplayFormat', '%.2e', 'Limits', [1e-14, 1e-1]);

            uilabel(gl, 'Text', 'Max step ("auto" or seconds):');
            if isnumeric(cfg.max_step) && ~isempty(cfg.max_step) && isfinite(cfg.max_step)
                msStr = num2str(cfg.max_step, '%.6g');
            else
                msStr = 'auto';
            end
            obj.edtMaxStep = uieditfield(gl, 'text', 'Value', msStr);

            uilabel(gl, 'Text', 'Secondary control:');
            scItems = {'No SC', 'AGC', 'Consensus', 'Passivity'};
            curSC = char(string(cfg.sc_structure));
            if ~ismember(curSC, scItems), scItems{end+1} = curSC; end
            obj.ddSC = uidropdown(gl, 'Items', scItems, 'Value', curSC);

            uilabel(gl, 'Text', 'Tertiary control:');
            tcItems = {'No TC', 'MAS'};
            curTC = char(string(cfg.tc_structure));
            if ~ismember(curTC, tcItems), tcItems{end+1} = curTC; end
            obj.ddTC = uidropdown(gl, 'Items', tcItems, 'Value', curTC);

            bar = uigridlayout(gl, [1, 3]);
            bar.Layout.Row = 8;
            bar.Layout.Column = [1, 2];
            bar.ColumnWidth = {'1x', 100, 90};
            bar.Padding = [0, 6, 0, 0];

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
            sess.config.sim.t_end = double(obj.edtTend.Value);
            sess.config.sim.method = char(obj.ddMethod.Value);
            sess.config.sim.rtol = double(obj.edtRtol.Value);
            sess.config.sim.atol = double(obj.edtAtol.Value);
            msTxt = strtrim(char(obj.edtMaxStep.Value));
            msNum = str2double(msTxt);
            if ~isnan(msNum) && msNum > 0
                sess.config.sim.max_step = msNum;
            else
                sess.config.sim.max_step = 'auto';
            end
            sess.config.sim.sc_structure = char(obj.ddSC.Value);
            sess.config.sim.tc_structure = char(obj.ddTC.Value);
            sess.log('Updated Simulation settings (t_end=%.2f s, solver=%s).', ...
                sess.config.sim.t_end, sess.config.sim.method);
            if ~isempty(obj.app.simTab)
                obj.app.simTab.syncControlsFromSession();
            end
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
