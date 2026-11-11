classdef ModalSettingsDialog < handle
    %MODALSETTINGSDIALOG Modal analysis sorting & frequency-filter dialog.

    properties
        app
        fig
        ddSort
        edtMinFreq
    end

    methods
        function obj = ModalSettingsDialog(app)
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            cfg = obj.app.session.config.modal;
            obj.fig = uifigure( ...
                'Name', 'Modal Analysis Settings', ...
                'Position', [310, 270, 400, 195], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [3, 2]);
            gl.RowHeight = {32, 32, '1x'};
            gl.ColumnWidth = {195, '1x'};
            gl.Padding = [14, 14, 14, 14];
            gl.RowSpacing = 8;

            uilabel(gl, 'Text', 'Sort modes by:');
            items = {'damping', 'frequency', 'real'};
            cur = lower(char(string(cfg.sort_by)));
            if ~ismember(cur, items), cur = 'damping'; end
            obj.ddSort = uidropdown(gl, 'Items', items, 'Value', cur);

            uilabel(gl, 'Text', 'Minimum frequency filter (Hz):');
            obj.edtMinFreq = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.min_freq), 'Limits', [0, 1e4]);

            bar = uigridlayout(gl, [1, 3]);
            bar.Layout.Row = 3;
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
            sess.config.modal.sort_by = char(obj.ddSort.Value);
            sess.config.modal.min_freq = double(obj.edtMinFreq.Value);
            sess.log('Updated Modal settings (sort_by=%s, min_freq=%.3g Hz).', ...
                sess.config.modal.sort_by, sess.config.modal.min_freq);
            if ~isempty(sess.linResult) && ~isempty(obj.app.modalTab)
                vflexp_modal(sess, false);
                obj.app.modalTab.refresh();
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
