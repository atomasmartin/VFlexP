classdef VFlexPApp < handle
    %VFLEXPAPP Main 6-tab programmatic MATLAB App for VFlexP.
    %   Provides Static Analysis, System Diagram, State Relevance,
    %   Complexity, Simulation, and Modal Analysis backed by VFlexPSession
    %   and modular Simulink libraries (zero binary .mlapp files).

    properties
        session VFlexPSession
        fig
        tabGroup

        lblCaseBadge
        lblQuickStatus
        lblLogStrip

        staticTab
        diagramTab
        relTab
        compTab
        simTab
        modalTab
    end

    methods
        function obj = VFlexPApp(casePath)
            obj.session = VFlexPSession();
            obj.buildUI();
            if nargin >= 1 && ~isempty(casePath)
                obj.loadCase(casePath);
            else
                obj.notifyCaseChanged();
            end
        end

        function buildUI(obj)
            obj.fig = uifigure( ...
                'Name', 'VFlexP — Power System Simulation & Analysis', ...
                'Tag', 'VFlexPMainFigure', ...
                'Position', [90, 65, 1340, 820], ...
                'Color', [0.95, 0.96, 0.98], ...
                'CloseRequestFcn', @(~,~) obj.delete());
            setappdata(obj.fig, 'VFlexPAppInstance', obj);

            rootGrid = uigridlayout(obj.fig, [3, 1]);
            rootGrid.RowHeight = {48, '1x', 28};
            rootGrid.Padding = [8, 8, 8, 6];
            rootGrid.RowSpacing = 6;

            % 1. Top global header bar
            hdr = uipanel(rootGrid, 'BackgroundColor', [0.11, 0.17, 0.26], 'BorderType', 'none');
            gHdr = uigridlayout(hdr, [1, 8]);
            gHdr.ColumnWidth = {145, 220, '1x', 110, 110, 115, 125, 80};
            gHdr.Padding = [12, 6, 12, 6];
            gHdr.ColumnSpacing = 8;

            uilabel(gHdr, ...
                'Text', 'VFlexP Studio', ...
                'FontSize', 16, ...
                'FontWeight', 'bold', ...
                'FontColor', [1, 1, 1]);

            obj.lblCaseBadge = uilabel(gHdr, ...
                'Text', 'Case: Untitled', ...
                'FontWeight', 'bold', ...
                'FontColor', [0.68, 0.85, 1.0]);

            obj.lblQuickStatus = uilabel(gHdr, ...
                'Text', 'Modular Simulink Library Engine Ready', ...
                'FontColor', [0.80, 0.85, 0.92]);

            uibutton(gHdr, 'Text', 'Load Case...', ...
                'ButtonPushedFcn', @(~,~) obj.onLoadCaseDialog());
            uibutton(gHdr, 'Text', 'Save Case...', ...
                'ButtonPushedFcn', @(~,~) obj.onSaveCaseDialog());
            uibutton(gHdr, 'Text', 'Variables ($)', ...
                'ButtonPushedFcn', @(~,~) VariablesDialog(obj));
            uibutton(gHdr, 'Text', 'Open Simulink', ...
                'ButtonPushedFcn', @(~,~) obj.openSimulink());
            uibutton(gHdr, 'Text', 'Logs', ...
                'ButtonPushedFcn', @(~,~) LogDialog(obj));

            % 2. Six main workspace tabs (aligned with VFlexPy; no POD, no Inertia)
            obj.tabGroup = uitabgroup(rootGrid);

            t1 = uitab(obj.tabGroup, 'Title', 'Static Analysis');
            obj.staticTab = StaticTab(t1, obj);

            t2 = uitab(obj.tabGroup, 'Title', 'System Diagram');
            obj.diagramTab = DiagramTab(t2, obj);

            t3 = uitab(obj.tabGroup, 'Title', 'State Relevance');
            obj.relTab = RelevanceTab(t3, obj);

            t4 = uitab(obj.tabGroup, 'Title', 'Complexity');
            obj.compTab = ComplexityTab(t4, obj);

            t5 = uitab(obj.tabGroup, 'Title', 'Simulation');
            obj.simTab = SimulationTab(t5, obj);

            t6 = uitab(obj.tabGroup, 'Title', 'Modal Analysis');
            obj.modalTab = ModalTab(t6, obj);

            % 3. Bottom live log strip
            botStrip = uigridlayout(rootGrid, [1, 2]);
            botStrip.ColumnWidth = {'1x', 95};
            botStrip.Padding = [6, 0, 6, 0];

            obj.lblLogStrip = uilabel(botStrip, ...
                'Text', 'Ready.', ...
                'FontName', 'Consolas', ...
                'FontSize', 11, ...
                'FontColor', [0.22, 0.28, 0.36]);
            uibutton(botStrip, 'Text', 'View Logs', ...
                'ButtonPushedFcn', @(~,~) LogDialog(obj));
        end

        function loadCase(obj, casePath)
            vflexp_load_case(obj.session, casePath);
            try
                vflexp_run_pf(obj.session);
            catch ME
                obj.session.log('Initial PF skipped/failed: %s', ME.message);
            end
            obj.notifyCaseChanged();
        end

        function saved = saveCase(obj, savePath)
            saved = vflexp_save_case(obj.session, savePath);
            obj.lblCaseBadge.Text = sprintf('Case: %s', obj.session.caseName);
            obj.refreshLogStrip();
        end

        function res = runPF(obj)
            res = vflexp_run_pf(obj.session);
            obj.notifyPFUpdated();
        end

        function res = runOPF(obj)
            res = vflexp_run_opf(obj.session);
            obj.notifyPFUpdated();
        end

        function notifyCaseChanged(obj)
            sess = obj.session;
            obj.lblCaseBadge.Text = sprintf('Case: %s', sess.caseName);
            obj.fig.Name = sprintf('VFlexP — %s', sess.caseName);
            obj.updateHeaderKPIs();

            obj.staticTab.refresh();
            obj.diagramTab.refresh();
            obj.relTab.refresh();
            obj.compTab.refresh();
            obj.simTab.refresh();
            obj.modalTab.refresh();
            obj.refreshLogStrip();
        end

        function notifyPFUpdated(obj)
            obj.updateHeaderKPIs();
            obj.staticTab.refresh();
            obj.diagramTab.refresh();
            obj.compTab.refresh();
            obj.refreshLogStrip();
        end

        function notifyVariablesChanged(obj)
            obj.staticTab.refresh();
            obj.diagramTab.refresh();
            obj.relTab.refresh();
            obj.compTab.refresh();
            obj.modalTab.refresh();
            obj.refreshLogStrip();
        end

        function notifyComplexityChanged(obj)
            obj.staticTab.refresh();
            obj.diagramTab.refresh();
            obj.compTab.refresh();
            obj.relTab.refresh();
            obj.modalTab.refresh();
            obj.refreshLogStrip();
        end

        function updateHeaderKPIs(obj)
            sess = obj.session;
            nB = height(sess.tables.Bus);
            nG = height(sess.tables.Gen);
            nL = height(sess.tables.Line);
            nIC = height(sess.tables.IntConv);
            nVars = numel(fieldnames(sess.variables));
            dyrTag = '';
            if ~isempty(sess.dyrData) && isstruct(sess.dyrData) && isfield(sess.dyrData, 'records')
                dyrTag = sprintf(' | .dyr (%d rec)', numel(sess.dyrData.records));
            end
            obj.lblQuickStatus.Text = sprintf( ...
                '%d Buses | %d Units | %d Lines | %d IntConvs | %d $vars%s', ...
                nB, nG, nL, nIC, nVars, dyrTag);
        end

        function refreshLogStrip(obj)
            if isempty(obj.lblLogStrip) || ~isvalid(obj.lblLogStrip)
                return;
            end
            if isempty(obj.session.logs)
                obj.lblLogStrip.Text = 'Ready.';
            else
                obj.lblLogStrip.Text = char(obj.session.logs(end));
            end
        end

        function openSimulink(obj)
            sess = obj.session;
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            mdlPath = fullfile(repoRoot, 'Simulink', 'VFlexP.slx');
            if ~bdIsLoaded('VFlexP')
                load_system(mdlPath);
            end
            if isempty(sess.initResult)
                try
                    if isempty(sess.pfResult) || ~sess.pfResult.converged
                        vflexp_run_pf(sess);
                    end
                    vflexp_build_and_init(sess);
                catch ME
                    sess.log('Opening VFlexP.slx with workspace defaults (%s).', ME.message);
                end
            end
            open_system('VFlexP');
            sess.log('Opened Simulink diagram VFlexP.slx (ModelWorkspace synced).');
            obj.refreshLogStrip();
        end

        function onLoadCaseDialog(obj)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            startDir = fullfile(repoRoot, 'bus_systems');
            if ~isfolder(startDir), startDir = repoRoot; end
            [file, path] = uigetfile( ...
                {'*.xlsx;*.xls', 'VFlexP Case Workbooks (*.xlsx, *.xls)'}, ...
                'Load VFlexP Network Case', startDir);
            if isequal(file, 0)
                return;
            end
            obj.loadCase(fullfile(path, file));
        end

        function onSaveCaseDialog(obj)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            startFile = fullfile(repoRoot, 'bus_systems', [obj.session.caseName '.xlsx']);
            [file, path] = uiputfile( ...
                {'*.xlsx', 'VFlexP Case Workbook (*.xlsx)'}, ...
                'Save VFlexP Network Case', startFile);
            if isequal(file, 0)
                return;
            end
            obj.saveCase(fullfile(path, file));
        end

        function delete(obj)
            if ~isempty(obj.fig) && isvalid(obj.fig)
                delete(obj.fig);
            end
        end
    end
end
