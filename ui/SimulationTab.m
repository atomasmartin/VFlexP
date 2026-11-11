classdef SimulationTab < handle
    %SIMULATIONTAB Time-domain simulation, disturbance events & 2x2 interactive scopes tab.

    properties
        app
        parent
        edtTend
        ddSolver
        lblStatus

        tblEvents
        selectedEventRow double = 0

        ddScope cell = cell(1, 4)
        axScope cell = cell(1, 4)
    end

    methods
        function obj = SimulationTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            obj.buildUI();
            obj.syncControlsFromSession();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [2, 1]);
            gl.RowHeight = {42, '1x'};
            gl.Padding = [10, 10, 10, 10];
            gl.RowSpacing = 8;

            % Toolbar
            bar = uigridlayout(gl, [1, 9]);
            bar.ColumnWidth = {65, 80, 52, 105, 105, 130, 95, 155, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uilabel(bar, 'Text', 't_end (s):', 'HorizontalAlignment', 'right');
            obj.edtTend = uieditfield(bar, 'numeric', ...
                'Value', 10.0, ...
                'Limits', [1e-3, 1e5], ...
                'ValueChangedFcn', @(~,~) obj.onSimToolbarEdited());

            uilabel(bar, 'Text', 'Solver:', 'HorizontalAlignment', 'right');
            obj.ddSolver = uidropdown(bar, ...
                'Items', {'ode15s', 'ode23tb', 'ode23s', 'ode23t', 'ode45'}, ...
                'Value', 'ode15s', ...
                'ValueChangedFcn', @(~,~) obj.onSimToolbarEdited());

            uibutton(bar, 'Text', 'Initialise', ...
                'ButtonPushedFcn', @(~,~) obj.onInitialise());

            uibutton(bar, 'Text', 'Play (Simulate)', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.runSimulation());

            uibutton(bar, 'Text', 'Settings...', ...
                'ButtonPushedFcn', @(~,~) DynamicSettingsDialog(obj.app));

            uibutton(bar, 'Text', 'Open Simulink Model', ...
                'BackgroundColor', [0.22, 0.28, 0.36], ...
                'FontColor', [1, 1, 1], ...
                'ButtonPushedFcn', @(~,~) obj.app.openSimulink());

            obj.lblStatus = uilabel(bar, ...
                'Text', 'Configure events and click "Play (Simulate)".', ...
                'FontColor', [0.25, 0.32, 0.40]);

            % Split body: Left = Disturbance Event Table, Right = 2x2 Scopes
            body = uigridlayout(gl, [1, 2]);
            body.ColumnWidth = {430, '1x'};
            body.Padding = [0, 0, 0, 0];
            body.ColumnSpacing = 10;

            pEv = uipanel(body, 'Title', 'Temporal Disturbances & Events (t_ev)');
            gEv = uigridlayout(pEv, [2, 1]);
            gEv.RowHeight = {'1x', 34};
            gEv.Padding = [6, 6, 6, 6];
            gEv.RowSpacing = 6;

            targetPop = {'AC Bus', 'DC Bus', 'AC Generator', 'DC Generator', 'Converter', 'AC Line'};
            chanPop = {'dP_load_MW', 'dQ_load_MVAr', 'dP_pu', 'dQ_pu', 'dI_pu', ...
                       'dP_ref', 'dV_ref', 'dQ_ref', 'df_ref', 'close'};
            obj.tblEvents = uitable(gEv, ...
                'ColumnName', {'On', 't_ev (s)', 'Target Type', 'Target ID', 'Channel', 'Magnitude', 'Duration (s)'}, ...
                'ColumnEditable', true(1, 7), ...
                'ColumnFormat', {'logical', 'numeric', targetPop, 'char', chanPop, 'numeric', 'numeric'}, ...
                'RowName', [], ...
                'CellEditCallback', @(~,~) obj.syncEventsToSession(), ...
                'CellSelectionCallback', @(~, evt) obj.onEventSelected(evt));

            evBtns = uigridlayout(gEv, [1, 4]);
            evBtns.ColumnWidth = {115, 110, 90, '1x'};
            evBtns.Padding = [0, 0, 0, 0];
            uibutton(evBtns, 'Text', '+ Add Event', 'ButtonPushedFcn', @(~,~) obj.onAddEvent());
            uibutton(evBtns, 'Text', '- Remove',    'ButtonPushedFcn', @(~,~) obj.onRemoveEvent());
            uibutton(evBtns, 'Text', 'Clear All',  'ButtonPushedFcn', @(~,~) obj.onClearEvents());
            uilabel(evBtns, 'Text', 'Inf = step', 'FontColor', [0.4, 0.45, 0.52], 'HorizontalAlignment', 'right');

            % 2x2 Scope grid
            gScopes = uigridlayout(body, [2, 2]);
            gScopes.Padding = [0, 0, 0, 0];
            gScopes.RowSpacing = 8;
            gScopes.ColumnSpacing = 8;

            scopeSignals = { ...
                'AC Bus Voltage Vm (pu)', ...
                'Frequency (Hz)', ...
                'Active Power P (MW)', ...
                'Reactive Power Q (MVAr)', ...
                'DC Bus Voltage Vdc (pu)', ...
                'AC Bus Angle (deg)', ...
                'Unit Terminal Voltage (pu)'};
            defaults = scopeSignals(1:4);

            for k = 1:4
                pS = uipanel(gScopes);
                gs = uigridlayout(pS, [2, 1]);
                gs.RowHeight = {28, '1x'};
                gs.Padding = [6, 4, 6, 4];
                gs.RowSpacing = 4;

                obj.ddScope{k} = uidropdown(gs, ...
                    'Items', scopeSignals, ...
                    'Value', defaults{k}, ...
                    'ValueChangedFcn', @(~,~) obj.renderScopes());
                obj.axScope{k} = uiaxes(gs);
                grid(obj.axScope{k}, 'on');
                xlabel(obj.axScope{k}, 'Time (s)');
            end
        end

        function syncControlsFromSession(obj)
            cfg = obj.app.session.config.sim;
            obj.edtTend.Value = double(cfg.t_end);
            sol = char(string(cfg.method));
            if ~ismember(sol, obj.ddSolver.Items)
                obj.ddSolver.Items = [obj.ddSolver.Items, {sol}];
            end
            obj.ddSolver.Value = sol;
            obj.populateEventsTable();
        end

        function refresh(obj)
            obj.syncControlsFromSession();
            obj.renderScopes();
        end

        function populateEventsTable(obj)
            ev = obj.app.session.eventList;
            if isempty(ev) || height(ev) == 0
                obj.tblEvents.Data = cell(0, 7);
                return;
            end
            n = height(ev);
            rows = cell(n, 7);
            for r = 1:n
                rows(r, :) = {logical(ev.Enabled(r)), double(ev.Time(r)), ...
                    char(ev.TargetType(r)), char(ev.TargetId(r)), ...
                    char(ev.Channel(r)), double(ev.Magnitude(r)), double(ev.Duration(r))};
            end
            obj.tblEvents.Data = rows;
        end

        function syncEventsToSession(obj)
            d = obj.tblEvents.Data;
            if isempty(d)
                obj.app.session.clearEvents();
                return;
            end
            n = size(d, 1);
            en = false(n, 1);
            tm = zeros(n, 1);
            tt = strings(n, 1);
            ti = strings(n, 1);
            ch = strings(n, 1);
            mg = zeros(n, 1);
            du = zeros(n, 1);
            for r = 1:n
                en(r) = logical(d{r, 1});
                tm(r) = double(d{r, 2});
                tt(r) = string(d{r, 3});
                ti(r) = string(d{r, 4});
                ch(r) = string(d{r, 5});
                mg(r) = double(d{r, 6});
                dv = double(d{r, 7});
                if isnan(dv) || dv <= 0
                    dv = Inf;
                end
                du(r) = dv;
            end
            obj.app.session.eventList = table(en, tm, tt, ti, ch, mg, du, ...
                'VariableNames', {'Enabled', 'Time', 'TargetType', 'TargetId', 'Channel', 'Magnitude', 'Duration'});
        end

        function onEventSelected(obj, evt)
            if ~isempty(evt.Indices)
                obj.selectedEventRow = evt.Indices(1, 1);
            end
        end

        function onAddEvent(obj)
            sess = obj.app.session;
            bId = '1';
            if ~isempty(sess.tables.Bus) && height(sess.tables.Bus) > 0
                bId = num2str(round(sess.tables.Bus.bus_i(1)));
            end
            tDefault = min(1.0, 0.2 * double(sess.config.sim.t_end));
            sess.addEvent(tDefault, 'AC Bus', bId, 'dP_load_MW', 5.0, Inf, true);
            obj.populateEventsTable();
            sess.log('Added disturbance event at t = %.2f s on AC Bus %s.', tDefault, bId);
            obj.app.refreshLogStrip();
        end

        function onRemoveEvent(obj)
            d = obj.tblEvents.Data;
            if isempty(d), return; end
            r = obj.selectedEventRow;
            if r < 1 || r > size(d, 1)
                r = size(d, 1);
            end
            d(r, :) = [];
            obj.tblEvents.Data = d;
            obj.syncEventsToSession();
        end

        function onClearEvents(obj)
            obj.app.session.clearEvents();
            obj.populateEventsTable();
        end

        function onSimToolbarEdited(obj)
            sess = obj.app.session;
            sess.config.sim.t_end = double(obj.edtTend.Value);
            sess.config.sim.method = char(obj.ddSolver.Value);
        end

        function initRes = onInitialise(obj)
            sess = obj.app.session;
            if isempty(sess.pfResult) || ~sess.pfResult.converged
                vflexp_run_pf(sess);
            end
            initRes = vflexp_build_and_init(sess);
            obj.lblStatus.Text = sprintf('Initialised (%d active states, %.2f s).', ...
                initRes.nStates, initRes.elapsed);
            obj.app.refreshLogStrip();
        end

        function simRes = runSimulation(obj)
            sess = obj.app.session;
            obj.onSimToolbarEdited();
            obj.syncEventsToSession();
            simRes = vflexp_simulate(sess, true);
            obj.lblStatus.Text = sprintf('Simulated [0, %.2f] s (%d steps) in %.2f s.', ...
                sess.config.sim.t_end, numel(simRes.t), simRes.elapsed);
            obj.renderScopes();
            obj.app.refreshLogStrip();
        end

        function renderScopes(obj)
            res = obj.app.session.simResult;
            for k = 1:4
                ax = obj.axScope{k};
                cla(ax);
                sigType = char(obj.ddScope{k}.Value);
                if isempty(res)
                    title(ax, sigType);
                    xlabel(ax, 'Time (s)');
                    grid(ax, 'on');
                    continue;
                end

                t = res.t;
                switch sigType
                    case 'AC Bus Voltage Vm (pu)'
                        Y = res.vac_mag;
                        lbls = res.busLabels;
                        yl = 'Vm (pu)';
                    case 'AC Bus Angle (deg)'
                        Y = res.vac_ang_deg;
                        lbls = res.busLabels;
                        yl = 'Angle (deg)';
                    case 'DC Bus Voltage Vdc (pu)'
                        Y = res.vdc;
                        lbls = res.dcBusLabels;
                        yl = 'Vdc (pu)';
                    case 'Frequency (Hz)'
                        Y = res.omega_Hz;
                        lbls = res.agentLabels;
                        yl = 'f (Hz)';
                    case 'Active Power P (MW)'
                        Y = res.P_MW;
                        lbls = res.agentLabels;
                        yl = 'P (MW)';
                    case 'Reactive Power Q (MVAr)'
                        Y = res.Q_MVAr;
                        lbls = res.agentLabels;
                        yl = 'Q (MVAr)';
                    otherwise % 'Unit Terminal Voltage (pu)'
                        Y = res.Vgen_pu;
                        lbls = res.agentLabels;
                        yl = 'V (pu)';
                end

                if isempty(Y)
                    title(ax, [sigType ' (none in case)']);
                else
                    plot(ax, t, Y, 'LineWidth', 1.4);
                    ylabel(ax, yl);
                    xlabel(ax, 'Time (s)');
                    grid(ax, 'on');
                    if numel(lbls) <= 8 && ~isempty(lbls)
                        legend(ax, cellstr(lbls), 'Location', 'best', 'Interpreter', 'none', 'FontSize', 8);
                    end
                    title(ax, sigType);
                end
            end
        end
    end
end
