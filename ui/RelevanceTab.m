classdef RelevanceTab < handle
    %RELEVANCETAB State Relevance analysis, model-reduction & parametric sweeps tab.
    %   Implements the state-relevance methodology of Tomas-Martin et al. (2023)
    %   with one-click complexity reduction ("Set complexity") and dynamic/static sweeps.

    properties
        app
        parent
        subTabs

        edtThreshold
        lblSummary
        axRel
        tblRank
        tblPlan

        sweepDynamicTab
        sweepStaticTab
    end

    methods
        function obj = RelevanceTab(parent, app)
            obj.parent = parent;
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            gl = uigridlayout(obj.parent, [1, 1]);
            gl.Padding = [6, 6, 6, 6];

            obj.subTabs = uitabgroup(gl);

            % 1. Single analysis sub-tab
            tSingle = uitab(obj.subTabs, 'Title', 'Single analysis');
            gS = uigridlayout(tSingle, [3, 1]);
            gS.RowHeight = {40, '1.15x', '1x'};
            gS.Padding = [8, 8, 8, 8];
            gS.RowSpacing = 8;

            bar = uigridlayout(gS, [1, 6]);
            bar.ColumnWidth = {150, 75, 95, 135, 100, '1x'};
            bar.Padding = [0, 0, 0, 0];
            bar.ColumnSpacing = 8;

            uibutton(bar, 'Text', 'Compute Relevance', ...
                'BackgroundColor', [0.13, 0.45, 0.82], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.runSingle());

            uilabel(bar, 'Text', 'Threshold:', 'HorizontalAlignment', 'right');
            obj.edtThreshold = uieditfield(bar, 'numeric', ...
                'Value', double(obj.app.session.config.relevance.threshold), ...
                'ValueDisplayFormat', '%.4g', ...
                'Limits', [0, 1], ...
                'ValueChangedFcn', @(~,~) obj.onThresholdChanged());

            uibutton(bar, 'Text', 'Set complexity', ...
                'BackgroundColor', [0.16, 0.56, 0.36], ...
                'FontColor', [1, 1, 1], ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.applySimplification());

            uibutton(bar, 'Text', 'Settings...', ...
                'ButtonPushedFcn', @(~,~) RelevanceSettingsDialog(obj.app));

            obj.lblSummary = uilabel(bar, ...
                'Text', 'Click "Compute Relevance" to rank states via balanced realization + participation factors.', ...
                'FontColor', [0.25, 0.32, 0.40]);

            % Log-scale relevance plot
            obj.axRel = uiaxes(gS);
            title(obj.axRel, 'Ranked State Relevance \rho_i (Log Scale)');
            ylabel(obj.axRel, 'State Relevance \rho_i');
            grid(obj.axRel, 'on');

            % Bottom tables: ranked states + candidate complexity reductions
            bot = uigridlayout(gS, [1, 2]);
            bot.ColumnWidth = {'1.25x', '1x'};
            bot.Padding = [0, 0, 0, 0];
            bot.ColumnSpacing = 8;

            pRank = uipanel(bot, 'Title', 'Ranked Dynamic States');
            gr = uigridlayout(pRank, [1, 1]); gr.Padding = [4, 4, 4, 4];
            obj.tblRank = uitable(gr, ...
                'ColumnName', {'Rank', 'State', 'Domain', 'Element', 'Sub-group', 'Relevance', 'Accum. (%)', '< Thresh'}, ...
                'ColumnEditable', false(1, 8), ...
                'RowName', []);

            pPlan = uipanel(bot, 'Title', 'Candidate Simplifications (for "Set complexity")');
            gp = uigridlayout(pPlan, [1, 1]); gp.Padding = [4, 4, 4, 4];
            obj.tblPlan = uitable(gp, ...
                'ColumnName', {'Kind', 'Element', 'Current', 'Proposed', 'Max Relevance'}, ...
                'ColumnEditable', false(1, 5), ...
                'RowName', []);

            % 2. Sweep: dynamic model (keeps PF operating point)
            tSweepDyn = uitab(obj.subTabs, 'Title', 'Sweep: dynamic model');
            obj.sweepDynamicTab = SweepTab(tSweepDyn, obj.app, 'relevance', 'keep');

            % 3. Sweep: static model (recomputes PF at every step)
            tSweepSta = uitab(obj.subTabs, 'Title', 'Sweep: static model');
            obj.sweepStaticTab = SweepTab(tSweepSta, obj.app, 'relevance', 'recompute');
        end

        function syncControlsFromSession(obj)
            obj.edtThreshold.Value = double(obj.app.session.config.relevance.threshold);
            if ~isempty(obj.app.session.relevanceResult)
                obj.onThresholdChanged();
            end
        end

        function refresh(obj)
            obj.edtThreshold.Value = double(obj.app.session.config.relevance.threshold);
            res = obj.app.session.relevanceResult;
            if isempty(res)
                cla(obj.axRel);
                title(obj.axRel, 'Ranked State Relevance \rho_i (Log Scale)');
                obj.tblRank.Data = cell(0, 8);
                obj.tblPlan.Data = cell(0, 5);
                obj.lblSummary.Text = 'Click "Compute Relevance" to rank states and preview model simplifications.';
            else
                obj.renderSingle(res);
            end
            obj.sweepDynamicTab.refresh();
            obj.sweepStaticTab.refresh();
        end

        function onThresholdChanged(obj)
            sess = obj.app.session;
            thr = double(obj.edtThreshold.Value);
            sess.config.relevance.threshold = thr;
            if ~isempty(sess.relevanceResult)
                res = sess.relevanceResult;
                res.threshold = thr;
                res.rankTable.BelowThreshold = res.rankTable.Relevance < thr;
                res.planTable = vflexp_apply_relevance(sess, thr, res.stateRelevance, res.stateMeta, true);
                sess.relevanceResult = res;
                obj.renderSingle(res);
            end
        end

        function res = runSingle(obj)
            sess = obj.app.session;
            sess.config.relevance.threshold = double(obj.edtThreshold.Value);
            res = vflexp_relevance(sess, true);
            obj.renderSingle(res);
            obj.app.refreshLogStrip();
        end

        function appRes = applySimplification(obj)
            sess = obj.app.session;
            thr = double(obj.edtThreshold.Value);
            sess.config.relevance.threshold = thr;
            if isempty(sess.relevanceResult)
                vflexp_relevance(sess, true);
            end
            appRes = vflexp_apply_relevance(sess, thr);
            obj.lblSummary.Text = sprintf( ...
                'Applied %d complexity reduction(s) at threshold %.4g. Click "Compute Relevance" to re-evaluate reduced model.', ...
                appRes.nSimplified, thr);
            obj.app.notifyComplexityChanged();
        end

        function renderSingle(obj, res)
            thr = double(obj.edtThreshold.Value);
            rt = res.rankTable;
            n = height(rt);

            cla(obj.axRel);
            if n > 0
                vals = max(double(rt.Relevance(:)), 1e-12);
                below = vals < thr;
                hold(obj.axRel, 'on');
                aboveIdx = find(~below);
                belowIdx = find(below);
                if ~isempty(aboveIdx)
                    bar(obj.axRel, aboveIdx, vals(aboveIdx), 0.72, ...
                        'FaceColor', [0.16, 0.48, 0.84], 'EdgeColor', 'none', 'DisplayName', '>= Threshold');
                end
                if ~isempty(belowIdx)
                    bar(obj.axRel, belowIdx, vals(belowIdx), 0.72, ...
                        'FaceColor', [0.72, 0.78, 0.85], 'EdgeColor', 'none', 'DisplayName', '< Threshold (reducible)');
                end
                yline(obj.axRel, max(thr, 1e-12), '--r', sprintf('Threshold = %.3g', thr), 'LineWidth', 1.4);
                set(obj.axRel, 'YScale', 'log');
                if n <= 35
                    obj.axRel.XTick = 1:n;
                    obj.axRel.XTickLabel = cellstr(rt.State);
                    obj.axRel.XTickLabelRotation = 35;
                else
                    obj.axRel.XTickMode = 'auto';
                    obj.axRel.XTickLabelMode = 'auto';
                    xlabel(obj.axRel, 'Rank Index');
                end
                obj.axRel.TickLabelInterpreter = 'none';
                xlim(obj.axRel, [0.3, n + 0.7]);
                hold(obj.axRel, 'off');
                grid(obj.axRel, 'on');
                title(obj.axRel, sprintf('Ranked State Relevance (%d states, Lyapunov shift = %.3g)', n, res.shift));
                ylabel(obj.axRel, 'Relevance \rho_i (log)');
            end

            rankRows = cell(n, 8);
            for r = 1:n
                rankRows(r, :) = {rt.Rank(r), char(rt.State(r)), char(rt.Domain(r)), ...
                    char(rt.ElementKind(r)), char(rt.SubGroup(r)), rt.Relevance(r), ...
                    100.0 * rt.Accumulated(r), logical(rt.Relevance(r) < thr)};
            end
            obj.tblRank.Data = rankRows;

            pt = res.planTable;
            np = height(pt);
            planRows = cell(np, 5);
            for r = 1:np
                planRows(r, :) = {char(pt.Kind(r)), char(pt.Element(r)), ...
                    char(pt.OldComplexity(r)), char(pt.NewComplexity(r)), pt.MaxRelevance(r)};
            end
            obj.tblPlan.Data = planRows;

            obj.lblSummary.Text = sprintf( ...
                '%d states total  |  %d states < %.3g  |  %d model simplification(s) ready.', ...
                n, sum(rt.Relevance < thr), thr, np);
        end
    end
end
