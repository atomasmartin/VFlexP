classdef SweepTab < handle
    %SWEEPTAB Reusable parametric sweep panel for Modal and State Relevance tabs.
    %   Supports kind = 'modal' | 'relevance' and pfMode = 'keep' | 'recompute'.

    properties
        app
        parent
        kind char = 'modal'
        pfMode char = 'keep'

        ddTargetType
        ddParam
        edtFrom
        edtTo
        spnSteps
        btnRun
        lblStatus

        axMain
        axSecondary
        tblSummary
    end

    methods
        function obj = SweepTab(parent, app, kind, pfMode)
            obj.parent = parent;
            obj.app = app;
            obj.kind = lower(char(string(kind)));
            obj.pfMode = lower(char(string(pfMode)));
            obj.buildUI();
            obj.syncTargets();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [2, 1]);
            gl.RowHeight = {42, '1x'};
            gl.Padding = [8, 8, 8, 8];
            gl.RowSpacing = 8;

            % Toolbar
            bar = uigridlayout(gl, [1, 11]);
            bar.ColumnWidth = {52, 145, 72, 175, 42, 78, 28, 78, 45, 70, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 6;

            uilabel(bar, 'Text', 'Target:', 'HorizontalAlignment', 'right');
            obj.ddTargetType = uidropdown(bar, ...
                'Items', {'$Variable', 'Unit Parameter', 'Load Scale (Pd,Qd)', 'Table Cell'}, ...
                'ItemsData', {'variable', 'unitParam', 'loadScale', 'tableCell'}, ...
                'Value', 'variable', ...
                'ValueChangedFcn', @(~,~) obj.onTargetTypeChanged());

            uilabel(bar, 'Text', 'Parameter:', 'HorizontalAlignment', 'right');
            obj.ddParam = uidropdown(bar, ...
                'Items', {'Pd_Qd'}, ...
                'Value', 'Pd_Qd', ...
                'Editable', 'on', ...
                'ValueChangedFcn', @(~,~) obj.onParamSelected());

            uilabel(bar, 'Text', 'From:', 'HorizontalAlignment', 'right');
            obj.edtFrom = uieditfield(bar, 'numeric', 'Value', 0.8, 'ValueDisplayFormat', '%.4g');

            uilabel(bar, 'Text', 'To:', 'HorizontalAlignment', 'right');
            obj.edtTo = uieditfield(bar, 'numeric', 'Value', 1.2, 'ValueDisplayFormat', '%.4g');

            uilabel(bar, 'Text', 'Steps:', 'HorizontalAlignment', 'right');
            obj.spnSteps = uispinner(bar, 'Value', 5, 'Limits', [2, 200], 'Step', 1);

            rightBar = uigridlayout(bar, [1, 2]);
            rightBar.ColumnWidth = {110, '1x'};
            rightBar.Padding = [0, 0, 0, 0];
            rightBar.ColumnSpacing = 8;

            obj.btnRun = uibutton(rightBar, ...
                'Text', 'Run Sweep', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.runSweep());
            if strcmp(obj.pfMode, 'recompute')
                modeNote = 'Recomputes PF at every step (static + dynamic operating point)';
            else
                modeNote = 'Keeps fixed PF operating point (dynamic parameter sensitivity)';
            end
            obj.lblStatus = uilabel(rightBar, ...
                'Text', modeNote, ...
                'FontColor', [0.28, 0.34, 0.42]);

            % Content split
            if strcmp(obj.kind, 'modal')
                body = uigridlayout(gl, [2, 2]);
                body.RowHeight = {'1x', '1x'};
                body.ColumnWidth = {'1.15x', '1x'};
                body.Padding = [0, 0, 0, 0];
                body.RowSpacing = 8;
                body.ColumnSpacing = 8;

                obj.axMain = uiaxes(body);
                obj.axMain.Layout.Row = [1, 2];
                obj.axMain.Layout.Column = 1;
                title(obj.axMain, 'Eigenvalue Root Locus');
                xlabel(obj.axMain, 'Real Part \sigma (1/s)');
                ylabel(obj.axMain, 'Imaginary Part \omega (rad/s)');
                grid(obj.axMain, 'on');

                obj.axSecondary = uiaxes(body);
                obj.axSecondary.Layout.Row = 1;
                obj.axSecondary.Layout.Column = 2;
                title(obj.axSecondary, 'Minimum Damping Ratio vs Parameter');
                xlabel(obj.axSecondary, 'Parameter Value');
                ylabel(obj.axSecondary, 'Damping (%)');
                grid(obj.axSecondary, 'on');

                obj.tblSummary = uitable(body, ...
                    'ColumnName', {'Step', 'Value', 'Converged', 'Max Real (1/s)', 'Min Damping (%)', 'Crit. Freq (Hz)'}, ...
                    'ColumnEditable', false(1, 6), ...
                    'RowName', []);
                obj.tblSummary.Layout.Row = 2;
                obj.tblSummary.Layout.Column = 2;
            else
                body = uigridlayout(gl, [1, 2]);
                body.ColumnWidth = {'1.25x', '1x'};
                body.Padding = [0, 0, 0, 0];
                body.ColumnSpacing = 8;

                obj.axMain = uiaxes(body);
                title(obj.axMain, 'State Relevance Trajectories vs Parameter');
                xlabel(obj.axMain, 'Parameter Value');
                ylabel(obj.axMain, 'State Relevance \rho_i');
                grid(obj.axMain, 'on');

                obj.tblSummary = uitable(body, ...
                    'ColumnName', {'Step', 'Value', 'Converged', 'States < Thresh', 'Lyap. Shift', 'Top State'}, ...
                    'ColumnEditable', false(1, 6), ...
                    'RowName', []);
            end
        end

        function refresh(obj)
            obj.syncTargets();
            res = obj.getCachedResult();
            if ~isempty(res)
                obj.renderResult(res);
            end
        end

        function syncTargets(obj)
            sess = obj.app.session;
            fn = fieldnames(sess.variables);
            if strcmp(obj.ddTargetType.Value, 'variable') && isempty(fn)
                obj.ddTargetType.Value = 'unitParam';
            elseif ~isempty(fn) && strcmp(obj.ddTargetType.Value, 'unitParam') && isempty(obj.getCachedResult())
                obj.ddTargetType.Value = 'variable';
            end
            obj.onTargetTypeChanged();
        end

        function onTargetTypeChanged(obj)
            sess = obj.app.session;
            tt = char(obj.ddTargetType.Value);
            switch tt
                case 'variable'
                    fn = fieldnames(sess.variables);
                    if isempty(fn)
                        items = {'var_1'};
                    else
                        items = fn(:)';
                    end
                    obj.ddParam.Items = items;
                    obj.ddParam.Value = items{1};
                    obj.onParamSelected();

                case 'unitParam'
                    cands = {};
                    if ~isempty(sess.tables.Gen) && height(sess.tables.Gen) > 0
                        uTypes = unique(cellstr(sess.tables.Gen.type), 'stable');
                        for k = 1:numel(uTypes)
                            ut = uTypes{k};
                            switch ut
                                case 'VSM'
                                    cands = [cands, {'VSM.Hv', 'VSM.Dv', 'VSM.VI_L'}]; %#ok<AGROW>
                                case 'GFM'
                                    cands = [cands, {'GFM.mp', 'GFM.nq', 'GFM.Tf', 'GFM.VI_L'}]; %#ok<AGROW>
                                case 'GFL'
                                    cands = [cands, {'GFL.Kp_PLL', 'GFL.Ki_PLL', 'GFL.Kp_CC'}]; %#ok<AGROW>
                                case {'SG', 'MGES'}
                                    cands = [cands, {[ut '.H'], [ut '.D'], [ut '.R_droop']}]; %#ok<AGROW>
                                case 'DC_PES'
                                    cands = [cands, {'DC_PES.R_droop', 'DC_PES.Tau'}]; %#ok<AGROW>
                            end
                        end
                    end
                    if ~isempty(sess.tables.IntConv) && height(sess.tables.IntConv) > 0
                        icTypes = unique(cellstr(sess.tables.IntConv.type), 'stable');
                        for k = 1:numel(icTypes)
                            ut = icTypes{k};
                            if contains(ut, 'GFL')
                                cands = [cands, {[ut '.Kp_PLL'], [ut '.Kp_dc']}]; %#ok<AGROW>
                            elseif contains(ut, 'VSM')
                                cands = [cands, {[ut '.Hv'], [ut '.Dv']}]; %#ok<AGROW>
                            elseif contains(ut, 'GFM')
                                cands = [cands, {[ut '.mp'], [ut '.nq']}]; %#ok<AGROW>
                            end
                        end
                    end
                    if isempty(cands)
                        cands = {'VSM.Hv', 'GFM.mp', 'GFL.Kp_PLL'};
                    end
                    cands = unique(cands, 'stable');
                    obj.ddParam.Items = cands;
                    obj.ddParam.Value = cands{1};
                    obj.onParamSelected();

                case 'loadScale'
                    obj.ddParam.Items = {'Pd_Qd'};
                    obj.ddParam.Value = 'Pd_Qd';
                    obj.edtFrom.Value = 0.8;
                    obj.edtTo.Value = 1.2;
                    obj.spnSteps.Value = 5;

                case 'tableCell'
                    cands = {'Gen.Pg(1)', 'Bus.Pd(1)', 'Line.x(1)'};
                    obj.ddParam.Items = cands;
                    obj.ddParam.Value = cands{1};
                    obj.onParamSelected();
            end
        end

        function onParamSelected(obj)
            sess = obj.app.session;
            tt = char(obj.ddTargetType.Value);
            pName = strtrim(char(obj.ddParam.Value));
            if strcmp(tt, 'variable')
                cleanP = regexprep(pName, '^\$', '');
                vr = sess.variableRanges;
                if ~isempty(vr) && height(vr) > 0
                    idx = find(strcmp(string(vr.Parameter), string(cleanP)), 1);
                    if ~isempty(idx)
                        obj.edtFrom.Value = double(vr.From(idx));
                        obj.edtTo.Value = double(vr.To(idx));
                        obj.spnSteps.Value = max(2, round(double(vr.N_steps(idx))));
                        return;
                    end
                end
                if isfield(sess.variables, cleanP)
                    v = double(sess.variables.(cleanP));
                    if v ~= 0
                        obj.edtFrom.Value = 0.8 * v;
                        obj.edtTo.Value = 1.2 * v;
                    else
                        obj.edtFrom.Value = 0.0;
                        obj.edtTo.Value = 1.0;
                    end
                end
            elseif strcmp(tt, 'unitParam')
                tok = regexp(pName, '^([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)$', 'tokens', 'once');
                if ~isempty(tok)
                    repoRoot = fileparts(fileparts(mfilename('fullpath')));
                    S = load(fullfile(repoRoot, 'Simulink', 'lib', 'vflexp_lib_manifest.mat'), 'manifest');
                    if isfield(S.manifest.units, tok{1}) && isfield(S.manifest.units.(tok{1}).defaultPrm, tok{2})
                        dv = S.manifest.units.(tok{1}).defaultPrm.(tok{2});
                        if isnumeric(dv) && isscalar(dv) && dv ~= 0
                            obj.edtFrom.Value = 0.5 * double(dv);
                            obj.edtTo.Value = 1.5 * double(dv);
                        end
                    end
                end
            elseif strcmp(tt, 'tableCell')
                tok = regexp(pName, '^([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\((\d+)\)$', 'tokens', 'once');
                if ~isempty(tok) && isfield(sess.tables, tok{1})
                    T = sess.tables.(tok{1});
                    rIdx = str2double(tok{3});
                    if ismember(tok{2}, T.Properties.VariableNames) && height(T) >= rIdx
                        v = double(T.(tok{2})(rIdx));
                        if v ~= 0
                            obj.edtFrom.Value = 0.8 * v;
                            obj.edtTo.Value = 1.2 * v;
                        end
                    end
                end
            end
        end

        function res = runSweep(obj)
            sess = obj.app.session;
            spec = struct( ...
                'kind', obj.kind, ...
                'mode', obj.pfMode, ...
                'targetType', char(obj.ddTargetType.Value), ...
                'paramName', strtrim(char(obj.ddParam.Value)), ...
                'from', double(obj.edtFrom.Value), ...
                'to', double(obj.edtTo.Value), ...
                'nSteps', round(double(obj.spnSteps.Value)));

            obj.btnRun.Enable = 'off';
            cleanupBtn = onCleanup(@() set(obj.btnRun, 'Enable', 'on'));
            if isempty(cleanupBtn), res = []; return; end

            try
                res = vflexp_sweep(sess, spec, @(s, n, v) obj.updateProgress(s, n, v));
                obj.renderResult(res);
                obj.lblStatus.Text = sprintf('Completed %d/%d steps in %.2f s (%s).', ...
                    sum(res.converged), numel(res.values), res.elapsed, res.paramName);
                obj.app.refreshLogStrip();
            catch ME
                obj.lblStatus.Text = ['Sweep error: ' ME.message];
                sess.log('ERROR in %s sweep: %s', obj.kind, ME.message);
                obj.app.refreshLogStrip();
                rethrow(ME);
            end
        end

        function updateProgress(obj, stepIdx, nSteps, val)
            if ~isempty(obj.lblStatus) && isvalid(obj.lblStatus)
                obj.lblStatus.Text = sprintf('Running step %d/%d (%s = %.4g)...', ...
                    stepIdx, nSteps, strtrim(char(obj.ddParam.Value)), val);
                drawnow limitrate;
            end
        end

        function renderResult(obj, res)
            if isempty(res)
                return;
            end
            vals = res.values(:);
            nSteps = numel(vals);

            if strcmp(obj.kind, 'modal')
                cla(obj.axMain);
                hold(obj.axMain, 'on');
                eigMat = res.eigenvalueTraj;
                cmap = parula(max(nSteps, 2));
                for s = 1:nSteps
                    if res.converged(s)
                        lam = eigMat(:, s);
                        mk = 'o';
                        ms = 28;
                        if s == 1
                            mk = 's'; ms = 42;
                        elseif s == nSteps
                            mk = '^'; ms = 46;
                        end
                        scatter(obj.axMain, real(lam), imag(lam), ms, repmat(cmap(s, :), numel(lam), 1), ...
                            'filled', 'Marker', mk, 'MarkerEdgeColor', [0.15, 0.15, 0.15]);
                    end
                end
                xline(obj.axMain, 0, '--r', 'LineWidth', 1.1);
                hold(obj.axMain, 'off');
                grid(obj.axMain, 'on');
                title(obj.axMain, sprintf('Root Locus — %s ([%.3g .. %.3g], square=start, triangle=end)', ...
                    res.paramName, vals(1), vals(end)), 'Interpreter', 'none');
                xlabel(obj.axMain, 'Real Part (1/s)');
                ylabel(obj.axMain, 'Imaginary Part (rad/s)');

                cla(obj.axSecondary);
                minDamp = nan(nSteps, 1);
                rows = cell(nSteps, 6);
                for s = 1:nSteps
                    if res.converged(s) && ~isempty(res.stepResults{s})
                        sr = res.stepResults{s};
                        [mnD, mnIdx] = min(sr.dampingPct);
                        mxR = max(sr.realPart);
                        crF = sr.freqHz(mnIdx);
                        minDamp(s) = mnD;
                        rows(s, :) = {s, vals(s), 'Yes', mxR, mnD, crF};
                    else
                        rows(s, :) = {s, vals(s), char(res.errors(s)), NaN, NaN, NaN};
                    end
                end
                plot(obj.axSecondary, vals, minDamp, '-o', ...
                    'LineWidth', 1.8, 'Color', [0.13, 0.45, 0.82], 'MarkerFaceColor', [0.13, 0.45, 0.82]);
                yline(obj.axSecondary, 5.0, '--r', '5% target');
                grid(obj.axSecondary, 'on');
                title(obj.axSecondary, sprintf('Minimum Damping vs %s', res.paramName), 'Interpreter', 'none');
                xlabel(obj.axSecondary, res.paramName, 'Interpreter', 'none');
                ylabel(obj.axSecondary, 'Min Damping (%)');
                obj.tblSummary.Data = rows;
            else
                cla(obj.axMain);
                hold(obj.axMain, 'on');
                relMat = max(res.relevanceMatrix, 1e-12);
                nShow = min(size(relMat, 1), 16);
                meanRel = mean(relMat, 2, 'omitnan');
                [~, topOrd] = sort(meanRel, 'descend');
                showIdx = topOrd(1:nShow);
                cols = lines(max(nShow, 1));
                for k = 1:numel(showIdx)
                    stIdx = showIdx(k);
                    plot(obj.axMain, vals, relMat(stIdx, :), '-o', ...
                        'LineWidth', 1.4, 'MarkerSize', 4, 'Color', cols(k, :), ...
                        'DisplayName', char(res.stateLabels(stIdx)));
                end
                thr = obj.app.session.config.relevance.threshold;
                yline(obj.axMain, max(thr, 1e-12), '--r', sprintf('Threshold (%.3g)', thr), 'LineWidth', 1.2);
                set(obj.axMain, 'YScale', 'log');
                hold(obj.axMain, 'off');
                grid(obj.axMain, 'on');
                if nShow > 0
                    legend(obj.axMain, 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 8);
                end
                title(obj.axMain, sprintf('Top State Relevance Curves vs %s', res.paramName), 'Interpreter', 'none');
                xlabel(obj.axMain, res.paramName, 'Interpreter', 'none');
                ylabel(obj.axMain, 'State Relevance (log scale)');

                rows = cell(nSteps, 6);
                for s = 1:nSteps
                    if res.converged(s) && ~isempty(res.stepResults{s})
                        sr = res.stepResults{s};
                        nBelow = sum(sr.stateRelevance < thr);
                        topSt = char(sr.rankTable.State(1));
                        rows(s, :) = {s, vals(s), 'Yes', sprintf('%d / %d', nBelow, numel(sr.stateRelevance)), ...
                            res.shifts(s), topSt};
                    else
                        rows(s, :) = {s, vals(s), char(res.errors(s)), '-', NaN, '-'};
                    end
                end
                obj.tblSummary.Data = rows;
            end
        end

        function res = getCachedResult(obj)
            sess = obj.app.session;
            if strcmp(obj.kind, 'modal')
                if strcmp(obj.pfMode, 'recompute')
                    res = sess.sweepModalStatic;
                else
                    res = sess.sweepModalDynamic;
                end
            else
                if strcmp(obj.pfMode, 'recompute')
                    res = sess.sweepRelStatic;
                else
                    res = sess.sweepRelDynamic;
                end
            end
        end
    end
end
