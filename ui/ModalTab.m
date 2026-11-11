classdef ModalTab < handle
    %MODALTAB Small-signal modal analysis, participation factors & sweeps tab.
    %   Provides Single linearisation + Sweep: dynamic model + Sweep: static model.

    properties
        app
        parent
        subTabs

        ddPlotView
        lblSummary
        axModes
        tblModes
        axPart
        tblPart
        selectedModeIdx double = 1

        sweepDynamicTab
        sweepStaticTab
    end

    methods
        function obj = ModalTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [1, 1]);
            gl.Padding = [6, 6, 6, 6];

            obj.subTabs = uitabgroup(gl);

            % 1. Single linearisation sub-tab
            tSingle = uitab(obj.subTabs, 'Title', 'Single linearisation');
            gS = uigridlayout(tSingle, [2, 1]);
            gS.RowHeight = {40, '1x'};
            gS.Padding = [8, 8, 8, 8];
            gS.RowSpacing = 8;

            bar = uigridlayout(gS, [1, 6]);
            bar.ColumnWidth = {145, 72, 185, 105, 110, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uibutton(bar, 'Text', 'Compute Modal', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.runSingle());

            uilabel(bar, 'Text', 'Plot view:', 'HorizontalAlignment', 'right');
            obj.ddPlotView = uidropdown(bar, ...
                'Items', {'Complex Plane (sigma, j*omega)', 'Damping (%) vs Frequency (Hz)'}, ...
                'Value', 'Complex Plane (sigma, j*omega)', ...
                'ValueChangedFcn', @(~,~) obj.renderPlot());

            uibutton(bar, 'Text', 'Settings...', ...
                'ButtonPushedFcn', @(~,~) ModalSettingsDialog(obj.app));

            uibutton(bar, 'Text', 'Export...', ...
                'ButtonPushedFcn', @(~,~) obj.onExportDialog());

            obj.lblSummary = uilabel(bar, ...
                'Text', 'Click "Compute Modal" to linearise VFlexP and inspect eigenvalues & participations.', ...
                'FontColor', [0.25, 0.32, 0.40]);

            % 2x2 body: Left = Eigenvalue Plot + Modes Table; Right = Participation Bar + Table
            body = uigridlayout(gS, [2, 2]);
            body.RowHeight = {'1.1x', '1x'};
            body.ColumnWidth = {'1.15x', '1x'};
            body.Padding = [0, 0, 0, 0];
            body.RowSpacing = 8;
            body.ColumnSpacing = 10;

            obj.axModes = uiaxes(body);
            obj.axModes.Layout.Row = 1;
            obj.axModes.Layout.Column = 1;
            grid(obj.axModes, 'on');
            title(obj.axModes, 'Eigenvalue Spectrum');

            obj.tblModes = uitable(body, ...
                'ColumnName', {'Mode', 'Real (1/s)', 'Imag (rad/s)', 'Freq (Hz)', 'Damping (%)', 'Dominant State', 'Dom. PF'}, ...
                'ColumnEditable', false(1, 7), ...
                'RowName', [], ...
                'CellSelectionCallback', @(~, evt) obj.onModeSelected(evt));
            obj.tblModes.Layout.Row = 2;
            obj.tblModes.Layout.Column = 1;

            obj.axPart = uiaxes(body);
            obj.axPart.Layout.Row = 1;
            obj.axPart.Layout.Column = 2;
            grid(obj.axPart, 'on');
            title(obj.axPart, 'Participation Factors (select a mode)');

            obj.tblPart = uitable(body, ...
                'ColumnName', {'State', 'Domain', 'Element', 'Participation Factor'}, ...
                'ColumnEditable', false(1, 4), ...
                'RowName', []);
            obj.tblPart.Layout.Row = 2;
            obj.tblPart.Layout.Column = 2;

            % 2. Sweep: dynamic model
            tSweepDyn = uitab(obj.subTabs, 'Title', 'Sweep: dynamic model');
            obj.sweepDynamicTab = SweepTab(tSweepDyn, obj.app, 'modal', 'keep');

            % 3. Sweep: static model
            tSweepSta = uitab(obj.subTabs, 'Title', 'Sweep: static model');
            obj.sweepStaticTab = SweepTab(tSweepSta, obj.app, 'modal', 'recompute');
        end

        function refresh(obj)
            res = obj.app.session.modalResult;
            if isempty(res)
                cla(obj.axModes);
                cla(obj.axPart);
                obj.tblModes.Data = cell(0, 7);
                obj.tblPart.Data = cell(0, 4);
                obj.lblSummary.Text = 'Click "Compute Modal" to linearise at the current operating point.';
            else
                obj.renderSingle(res);
            end
            obj.sweepDynamicTab.refresh();
            obj.sweepStaticTab.refresh();
        end

        function res = runSingle(obj)
            sess = obj.app.session;
            res = vflexp_modal(sess, true);
            obj.selectedModeIdx = 1;
            obj.renderSingle(res);
            obj.app.refreshLogStrip();
        end

        function renderSingle(obj, res)
            mt = res.modesTable;
            n = height(mt);
            rows = cell(n, 7);
            for r = 1:n
                rows(r, :) = {mt.Mode(r), mt.Real_1_s(r), mt.Imag_rad_s(r), ...
                    mt.Freq_Hz(r), mt.Damping_pct(r), char(mt.DominantState(r)), mt.DominantPF(r)};
            end
            obj.tblModes.Data = rows;

            minDamp = min(res.dampingPct);
            maxReal = max(res.realPart);
            nUnstable = sum(res.realPart > 1e-6);
            obj.lblSummary.Text = sprintf( ...
                '%d modes  |  Unstable: %d  |  Max Re = %.4g 1/s  |  Min damping = %.2f%% (%.3f s)', ...
                n, nUnstable, maxReal, minDamp, res.elapsed);

            if obj.selectedModeIdx < 1 || obj.selectedModeIdx > n
                obj.selectedModeIdx = 1;
            end
            obj.renderPlot();
            obj.renderParticipation(obj.selectedModeIdx);
        end

        function renderPlot(obj)
            res = obj.app.session.modalResult;
            cla(obj.axModes);
            if isempty(res)
                return;
            end
            viewMode = char(obj.ddPlotView.Value);
            hold(obj.axModes, 'on');
            if startsWith(viewMode, 'Complex')
                sig = res.realPart;
                omg = res.imagPart;
                isUnst = sig > 1e-6;
                if any(~isUnst)
                    scatter(obj.axModes, sig(~isUnst), omg(~isUnst), 48, [0.13, 0.45, 0.82], ...
                        'filled', 'DisplayName', 'Stable');
                end
                if any(isUnst)
                    scatter(obj.axModes, sig(isUnst), omg(isUnst), 64, [0.82, 0.16, 0.16], ...
                        'filled', 'Marker', 'd', 'DisplayName', 'Unstable');
                end
                if obj.selectedModeIdx >= 1 && obj.selectedModeIdx <= numel(sig)
                    m = obj.selectedModeIdx;
                    plot(obj.axModes, sig(m), omg(m), 'o', 'MarkerSize', 12, ...
                        'LineWidth', 2.0, 'Color', [0.92, 0.55, 0.08]);
                end
                xline(obj.axModes, 0, '--r', 'LineWidth', 1.2);
                % 5% damping rays
                yMax = max([abs(omg); 1]);
                zeta = 0.05;
                xRay = -yMax * (zeta / sqrt(1 - zeta^2));
                plot(obj.axModes, [0, xRay], [0,  yMax], ':k', 'HandleVisibility', 'off');
                plot(obj.axModes, [0, xRay], [0, -yMax], ':k', 'HandleVisibility', 'off');
                xlabel(obj.axModes, 'Real Part \sigma (1/s)');
                ylabel(obj.axModes, 'Imaginary Part \omega (rad/s)');
                title(obj.axModes, 'Eigenvalues in Complex Plane (dotted = 5% damping)');
            else
                scatter(obj.axModes, res.freqHz, res.dampingPct, 50, [0.13, 0.45, 0.82], 'filled');
                if obj.selectedModeIdx >= 1 && obj.selectedModeIdx <= numel(res.freqHz)
                    m = obj.selectedModeIdx;
                    plot(obj.axModes, res.freqHz(m), res.dampingPct(m), 'o', ...
                        'MarkerSize', 12, 'LineWidth', 2.0, 'Color', [0.92, 0.55, 0.08]);
                end
                yline(obj.axModes, 5.0, '--r', '5% damping');
                yline(obj.axModes, 0.0, '-r');
                xlabel(obj.axModes, 'Frequency (Hz)');
                ylabel(obj.axModes, 'Damping Ratio (%)');
                title(obj.axModes, 'Damping Ratio vs Oscillation Frequency');
            end
            hold(obj.axModes, 'off');
            grid(obj.axModes, 'on');
        end

        function onModeSelected(obj, evt)
            if isempty(evt.Indices)
                return;
            end
            obj.selectedModeIdx = evt.Indices(1, 1);
            obj.renderPlot();
            obj.renderParticipation(obj.selectedModeIdx);
        end

        function selectMode(obj, modeIdx)
            obj.selectedModeIdx = modeIdx;
            obj.renderPlot();
            obj.renderParticipation(modeIdx);
        end

        function renderParticipation(obj, modeIdx)
            res = obj.app.session.modalResult;
            cla(obj.axPart);
            if isempty(res) || modeIdx < 1 || modeIdx > size(res.participation, 2)
                return;
            end

            pfCol = res.participation(:, modeIdx);
            [sortedPF, ord] = sort(pfCol, 'descend');
            nShow = min(12, numel(sortedPF));
            topIdx = ord(1:nShow);

            barh(obj.axPart, flipud(sortedPF(1:nShow)), 'FaceColor', [0.16, 0.52, 0.82]);
            obj.axPart.YTick = 1:nShow;
            obj.axPart.YTickLabel = cellstr(flipud(res.stateLabels(topIdx)));
            obj.axPart.TickLabelInterpreter = 'none';
            xlabel(obj.axPart, 'Normalized Participation Factor');
            xlim(obj.axPart, [0, max(1.05 * max(sortedPF), 0.1)]);
            grid(obj.axPart, 'on');
            title(obj.axPart, sprintf('Mode %d (%.3f Hz, %.2f%% damping)', ...
                modeIdx, res.freqHz(modeIdx), res.dampingPct(modeIdx)));

            nAll = numel(ord);
            pRows = cell(nAll, 4);
            for k = 1:nAll
                idx = ord(k);
                pRows(k, :) = {char(res.stateLabels(idx)), char(res.stateMeta.Domain(idx)), ...
                    char(res.stateMeta.ElementKind(idx)), pfCol(idx)};
            end
            obj.tblPart.Data = pRows;
        end

        function exportModal(obj, outPath)
            sess = obj.app.session;
            vflexp_modal(sess, isempty(sess.modalResult), outPath);
            obj.app.refreshLogStrip();
        end

        function onExportDialog(obj)
            [file, path] = uiputfile( ...
                {'*.xlsx', 'Excel Workbook (*.xlsx)'; '*.mat', 'MATLAB MAT-file (*.mat)'}, ...
                'Export Modal Analysis Results', ...
                sprintf('%s_modal.xlsx', obj.app.session.caseName));
            if isequal(file, 0)
                return;
            end
            obj.exportModal(fullfile(path, file));
        end
    end
end
