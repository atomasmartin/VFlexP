classdef LogDialog < handle
    %LOGDIALOG Session log viewer window for VFlexPApp.

    properties
        app
        fig
        txtArea
    end

    methods
        function obj = LogDialog(app)
            obj.app = app;
            obj.buildUI();
            obj.refresh();
        end

        function buildUI(obj)
            obj.fig = uifigure( ...
                'Name', 'VFlexP — Session Logs', ...
                'Position', [220, 180, 720, 440], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [2, 1]);
            gl.RowHeight = {'1x', 38};
            gl.Padding = [12, 12, 12, 12];
            gl.RowSpacing = 8;

            obj.txtArea = uitextarea(gl, ...
                'Editable', 'off', ...
                'FontName', 'Consolas', ...
                'FontSize', 12);

            bar = uigridlayout(gl, [1, 4]);
            bar.ColumnWidth = {110, 110, '1x', 100};
            bar.Padding = [0, 0, 0, 0];

            uibutton(bar, 'Text', 'Refresh', ...
                'ButtonPushedFcn', @(~,~) obj.refresh());
            uibutton(bar, 'Text', 'Clear Logs', ...
                'ButtonPushedFcn', @(~,~) obj.onClear());
            uilabel(bar, 'Text', '');
            uibutton(bar, 'Text', 'Close', ...
                'ButtonPushedFcn', @(~,~) obj.close());
        end

        function refresh(obj)
            if isempty(obj.fig) || ~isvalid(obj.fig)
                return;
            end
            logs = obj.app.session.logs;
            if isempty(logs)
                obj.txtArea.Value = {'(No log entries yet)'};
            else
                obj.txtArea.Value = cellstr(logs);
            end
        end

        function onClear(obj)
            obj.app.session.clearLogs();
            obj.refresh();
            obj.app.refreshLogStrip();
        end

        function close(obj)
            if ~isempty(obj.fig) && isvalid(obj.fig)
                delete(obj.fig);
            end
        end
    end
end
