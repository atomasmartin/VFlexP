classdef RelevanceSettingsDialog < handle
    %RELEVANCESETTINGSDIALOG State-relevance I/O channel & Lyapunov margin dialog.

    properties
        app
        fig
        edtMargin
        edtThresh
        chkInBusP
        chkInBusQ
        chkInDcI
        chkInGenP
        chkInGenV
        chkOutVac
        chkOutVdc
        chkOutElec
    end

    methods
        function obj = RelevanceSettingsDialog(app)
            obj.app = app;
            obj.buildUI();
        end

        function buildUI(obj)
            cfg = obj.app.session.config.relevance;
            inC = cellstr(string(cfg.input_channels));
            outC = cellstr(string(cfg.output_channels));

            obj.fig = uifigure( ...
                'Name', 'State Relevance Settings', ...
                'Position', [280, 210, 480, 340], ...
                'Color', [0.97, 0.98, 0.99]);

            gl = uigridlayout(obj.fig, [5, 2]);
            gl.RowHeight = {30, 30, 88, 72, '1x'};
            gl.ColumnWidth = {220, '1x'};
            gl.Padding = [14, 14, 14, 14];
            gl.RowSpacing = 8;

            uilabel(gl, 'Text', 'Simplification threshold:');
            obj.edtThresh = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.threshold), ...
                'ValueDisplayFormat', '%.4g', ...
                'Limits', [0, 1]);

            uilabel(gl, 'Text', 'Lyapunov stability shift margin:');
            obj.edtMargin = uieditfield(gl, 'numeric', ...
                'Value', double(cfg.stability_margin), ...
                'ValueDisplayFormat', '%.2e', ...
                'Limits', [1e-6, 10]);

            pIn = uipanel(gl, 'Title', 'Disturbance Input Channels (B)');
            pIn.Layout.Row = 3;
            pIn.Layout.Column = [1, 2];
            gIn = uigridlayout(pIn, [2, 3]);
            gIn.Padding = [8, 4, 8, 4];
            obj.chkInBusP = uicheckbox(gIn, 'Text', 'AC Bus dP', 'Value', ismember('bus_P', inC));
            obj.chkInBusQ = uicheckbox(gIn, 'Text', 'AC Bus dQ', 'Value', ismember('bus_Q', inC));
            obj.chkInDcI  = uicheckbox(gIn, 'Text', 'DC Bus dI', 'Value', ismember('dc_I', inC));
            obj.chkInGenP = uicheckbox(gIn, 'Text', 'Gen/Conv dP*', 'Value', ismember('gen_P', inC));
            obj.chkInGenV = uicheckbox(gIn, 'Text', 'Gen/Conv dV*', 'Value', ismember('gen_V', inC));

            pOut = uipanel(gl, 'Title', 'Observed Output Channels (C)');
            pOut.Layout.Row = 4;
            pOut.Layout.Column = [1, 2];
            gOut = uigridlayout(pOut, [1, 3]);
            gOut.Padding = [8, 4, 8, 4];
            obj.chkOutVac  = uicheckbox(gOut, 'Text', 'AC Voltages (vd,vq)', 'Value', ismember('vac', outC));
            obj.chkOutVdc  = uicheckbox(gOut, 'Text', 'DC Voltages (Vdc)',   'Value', ismember('vdc', outC));
            obj.chkOutElec = uicheckbox(gOut, 'Text', 'Unit Telemetry (f,P,V,Q)', 'Value', ismember('elec', outC));

            bar = uigridlayout(gl, [1, 3]);
            bar.Layout.Row = 5;
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
            sess.config.relevance.threshold = double(obj.edtThresh.Value);
            sess.config.relevance.stability_margin = double(obj.edtMargin.Value);

            inC = {};
            if obj.chkInBusP.Value, inC{end+1} = 'bus_P'; end
            if obj.chkInBusQ.Value, inC{end+1} = 'bus_Q'; end
            if obj.chkInDcI.Value,  inC{end+1} = 'dc_I';  end
            if obj.chkInGenP.Value, inC{end+1} = 'gen_P'; end
            if obj.chkInGenV.Value, inC{end+1} = 'gen_V'; end
            if isempty(inC), inC = {'bus_P', 'gen_P'}; end
            sess.config.relevance.input_channels = inC;

            outC = {};
            if obj.chkOutVac.Value,  outC{end+1} = 'vac';  end
            if obj.chkOutVdc.Value,  outC{end+1} = 'vdc';  end
            if obj.chkOutElec.Value, outC{end+1} = 'elec'; end
            if isempty(outC), outC = {'vac', 'elec'}; end
            sess.config.relevance.output_channels = outC;

            sess.log('Updated State Relevance settings (threshold=%.4g, margin=%.2e).', ...
                sess.config.relevance.threshold, sess.config.relevance.stability_margin);
            if ~isempty(obj.app.relTab)
                obj.app.relTab.syncControlsFromSession();
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
