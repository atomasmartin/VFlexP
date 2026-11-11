classdef VFlexPSession < handle
    %VFLEXPSESSION Central in-memory state container for VFlexP.
    %   Holds network tables, $variables, .dyr dynamic records, .json
    %   settings, disturbance events, and results of PF/OPF, initialization,
    %   simulation, linearization, modal analysis, state relevance, and sweeps.
    %   Zero dependency on gcb, bdroot, or Simulink mask dialog tables.

    properties
        casePath char = ''
        caseName char = 'Untitled'
        rawTables struct = struct()
        tables struct = struct()
        variables struct = struct()
        variableRanges table = table()
        dyrData = []
        unitOverrides struct = struct()
        config struct = struct()
        eventList table = table()

        pfResult = []
        initResult = []
        simResult = []
        linResult = []
        modalResult = []
        relevanceResult = []

        sweepModalDynamic = []
        sweepModalStatic = []
        sweepRelDynamic = []
        sweepRelStatic = []

        logs string = strings(0, 1)
    end

    methods
        function obj = VFlexPSession(casePath, varargin)
            %VFLEXPSESSION Construct an empty or loaded VFlexP session.
            obj.resetDefaults();
            if nargin >= 1 && ~isempty(casePath)
                vflexp_load_case(obj, casePath, varargin{:});
            end
        end

        function resetDefaults(obj)
            %RESETDEFAULTS Initialize clean default tables, config, and events.
            obj.rawTables = VFlexPSession.emptyRawTables();
            obj.tables = VFlexPSession.emptyTables();
            obj.variables = struct();
            obj.variableRanges = table(strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
                'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
            obj.dyrData = [];
            obj.unitOverrides = struct();

            obj.config = struct( ...
                'model', struct( ...
                    'frame_policy', 'first-gen'), ...
                'powerflow', struct( ...
                    'tolerance', 1e-8, ...
                    'max_iterations', 50, ...
                    'enforce_q_limits', false, ...
                    'write_branch', false), ...
                'sim', struct( ...
                    't_end', 10.0, ...
                    'method', 'ode15s', ...
                    'rtol', 1e-4, ...
                    'atol', 1e-6, ...
                    'max_step', 'auto', ...
                    'sc_structure', 'No SC', ...
                    'define_sc_in_matlab', false, ...
                    'tc_structure', 'No TC'), ...
                'modal', struct( ...
                    'input_channels', {{'bus_P', 'bus_Q', 'gen_P', 'gen_V'}}, ...
                    'output_channels', {{'vac', 'vdc', 'elec'}}, ...
                    'sort_by', 'damping', ...
                    'min_freq', 0.0), ...
                'relevance', struct( ...
                    'input_channels', {{'bus_P', 'bus_Q', 'gen_P', 'gen_V'}}, ...
                    'output_channels', {{'vac', 'vdc', 'elec'}}, ...
                    'stability_margin', 1e-2, ...
                    'threshold', 0.02));

            obj.eventList = VFlexPSession.emptyEventsTable();
            obj.invalidateDynamicResults();
            obj.pfResult = [];
        end

        function invalidateDynamicResults(obj)
            %INVALIDATEDYNAMICRESULTS Clear cached dynamic/modal/relevance outputs.
            obj.initResult = [];
            obj.simResult = [];
            obj.linResult = [];
            obj.modalResult = [];
            obj.relevanceResult = [];
        end

        function log(obj, fmt, varargin)
            %LOG Append a timestamped entry to the session log.
            if nargin > 2
                msg = sprintf(fmt, varargin{:});
            else
                msg = char(string(fmt));
            end
            ts = char(datetime('now', 'Format', 'HH:mm:ss'));
            entry = string(sprintf('[%s] %s', ts, msg));
            obj.logs(end+1, 1) = entry;
            if numel(obj.logs) > 500
                obj.logs = obj.logs(end-499:end, 1);
            end
        end

        function clearLogs(obj)
            %CLEARLOGS Empty the session log buffer.
            obj.logs = strings(0, 1);
        end

        function syncVariablesToBase(obj)
            %SYNCVARIABLESTOBASE Mirror session $variables into MATLAB base workspace.
            fnames = fieldnames(obj.variables);
            for k = 1:numel(fnames)
                vn = fnames{k};
                if isvarname(vn)
                    assignin('base', vn, obj.variables.(vn));
                end
            end
        end

        function setVariable(obj, name, val)
            %SETVARIABLE Update a $variable value and refresh evaluated tables & .dyr.
            cleanName = regexprep(strtrim(char(string(name))), '^\$', '');
            if ~isvarname(cleanName)
                error('VFlexP:Session:InvalidVarName', 'Invalid variable name: "%s".', cleanName);
            end
            valNum = double(val);
            obj.variables.(cleanName) = valNum;
            assignin('base', cleanName, valNum);

            if ~isempty(obj.variableRanges) && ismember('Parameter', obj.variableRanges.Properties.VariableNames)
                pcol = regexprep(strtrim(string(obj.variableRanges.Parameter)), '^\$', '');
                idx = find(strcmp(pcol, cleanName), 1);
                if ~isempty(idx)
                    obj.variableRanges.Value(idx) = valNum;
                else
                    fromV = 0.8 * valNum;
                    toV = 1.2 * valNum;
                    if valNum == 0
                        fromV = 0.0;
                        toV = 1.0;
                    end
                    newRow = table(string(cleanName), valNum, fromV, toV, 5, ...
                        'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
                    obj.variableRanges = [obj.variableRanges; newRow];
                end
            end

            if ~isempty(obj.dyrData) && isstruct(obj.dyrData)
                if isfield(obj.dyrData, 'definitions') && isfield(obj.dyrData.definitions, cleanName)
                    obj.dyrData.definitions.(cleanName) = valNum;
                end
                if isfield(obj.dyrData, 'records')
                    for r = 1:numel(obj.dyrData.records)
                        refs = obj.dyrData.records(r).references;
                        if isstruct(refs)
                            rf = fieldnames(refs);
                            for j = 1:numel(rf)
                                if strcmp(refs.(rf{j}), cleanName)
                                    obj.dyrData.records(r).params.(rf{j}) = valNum;
                                end
                            end
                        end
                    end
                end
            end

            obj.evalTables();
            obj.invalidateDynamicResults();
        end

        function evalTables(obj)
            %EVALTABLES Evaluate rawTables ($variables and expressions) into numeric tables.
            obj.syncVariablesToBase();

            % 1. System base Sb and Fb table
            obj.tables.Sb = obj.evalScalar(obj.rawTables.Sb, 100);
            if isempty(obj.rawTables.Fb)
                obj.tables.Fb = table(1, 50, 'VariableNames', {'Area', 'Fb'});
            else
                RFb = obj.rawTables.Fb;
                nFb = height(RFb);
                areaCol = zeros(nFb, 1);
                fbCol = zeros(nFb, 1);
                for r = 1:nFb
                    areaCol(r) = obj.evalScalar(RFb.Area(r), 1);
                    fbCol(r) = obj.evalScalar(RFb.Fb(r), 50);
                end
                obj.tables.Fb = table(areaCol, fbCol, 'VariableNames', {'Area', 'Fb'});
            end

            % 2. Bus table
            RBus = obj.rawTables.Bus;
            if isempty(RBus)
                obj.tables.Bus = VFlexPSession.emptyTables().Bus;
            else
                numCols = {'bus_i','type','Pd','Qd','Gs','Bs','area','Vm','Va','baseKV','zone','Vmax','Vmin'};
                TBus = table();
                nB = height(RBus);
                for c = 1:numel(numCols)
                    cn = numCols{c};
                    if ismember(cn, RBus.Properties.VariableNames)
                        TBus.(cn) = obj.evalColumn(RBus.(cn), 0);
                    else
                        TBus.(cn) = zeros(nB, 1);
                    end
                end
                if ismember('Complexity_bus', RBus.Properties.VariableNames)
                    cb = strtrim(string(RBus.Complexity_bus));
                    cb(cb == "" | ismissing(cb)) = "static omega";
                    TBus.Complexity_bus = cb;
                else
                    TBus.Complexity_bus = repmat("static omega", nB, 1);
                end
                if ismember('Complexity_load', RBus.Properties.VariableNames)
                    cl = strtrim(string(RBus.Complexity_load));
                    for r = 1:nB
                        if cl(r) == "" || ismissing(cl(r))
                            if startsWith(TBus.Complexity_bus(r), "DC", 'IgnoreCase', true)
                                cl(r) = "DC ConstP";
                            else
                                cl(r) = "series RL static omega";
                            end
                        end
                    end
                    TBus.Complexity_load = cl;
                else
                    TBus.Complexity_load = repmat("series RL static omega", nB, 1);
                end
                for gc = {'Latitude', 'Longitude'}
                    gn = gc{1};
                    if ismember(gn, RBus.Properties.VariableNames)
                        TBus.(gn) = obj.evalColumn(RBus.(gn), 0);
                    else
                        TBus.(gn) = zeros(nB, 1);
                    end
                end
                obj.tables.Bus = TBus;
            end

            % 3. Line table
            RLine = obj.rawTables.Line;
            if isempty(RLine) || height(RLine) == 0
                obj.tables.Line = VFlexPSession.emptyTables().Line;
            else
                lineNumCols = {'fbus','tbus','r','x','b','rateA','rateB','rateC','ratio','angle','status','angmin','angmax'};
                TLine = table();
                nL = height(RLine);
                for c = 1:numel(lineNumCols)
                    cn = lineNumCols{c};
                    defV = 0;
                    if strcmp(cn, 'status'), defV = 1; end
                    if ismember(cn, RLine.Properties.VariableNames)
                        TLine.(cn) = obj.evalColumn(RLine.(cn), defV);
                    else
                        TLine.(cn) = repmat(defV, nL, 1);
                    end
                end
                if ismember('Complexity', RLine.Properties.VariableNames)
                    lc = strtrim(string(RLine.Complexity));
                    lc(lc == "" | ismissing(lc)) = "static omega";
                    TLine.Complexity = lc;
                else
                    TLine.Complexity = repmat("static omega", nL, 1);
                end
                obj.tables.Line = TLine;
            end

            % 4. Gen table
            RGen = obj.rawTables.Gen;
            if isempty(RGen) || height(RGen) == 0
                obj.tables.Gen = VFlexPSession.emptyTables().Gen;
            else
                genCols = VFlexPSession.genColumnNames();
                TGen = table();
                nG = height(RGen);
                for c = 1:numel(genCols)
                    cn = genCols{c};
                    if strcmp(cn, 'type')
                        if ismember('type', RGen.Properties.VariableNames)
                            tp = strtrim(string(RGen.type));
                            tp(tp == "" | ismissing(tp)) = "GFM";
                        else
                            tp = repmat("GFM", nG, 1);
                        end
                        TGen.type = VFlexPSession.normalizeGenTypes(tp);
                    else
                        defV = 0;
                        if strcmp(cn, 'Vg') || strcmp(cn, 'status'), defV = 1; end
                        if strcmp(cn, 'mBase'), defV = obj.tables.Sb; end
                        if ismember(cn, RGen.Properties.VariableNames)
                            TGen.(cn) = obj.evalColumn(RGen.(cn), defV);
                        else
                            TGen.(cn) = repmat(defV, nG, 1);
                        end
                    end
                end
                obj.tables.Gen = TGen;
            end

            % 5. IntConv table
            RIC = obj.rawTables.IntConv;
            if isempty(RIC) || height(RIC) == 0
                obj.tables.IntConv = VFlexPSession.emptyTables().IntConv;
            else
                icCols = VFlexPSession.intConvColumnNames();
                TIC = table();
                nIC = height(RIC);
                for c = 1:numel(icCols)
                    cn = icCols{c};
                    if strcmp(cn, 'type')
                        if ismember('type', RIC.Properties.VariableNames)
                            tp = strtrim(string(RIC.type));
                            tp(tp == "" | ismissing(tp)) = "intconv_GFL";
                        else
                            tp = repmat("intconv_GFL", nIC, 1);
                        end
                        TIC.type = VFlexPSession.normalizeIntConvTypes(tp);
                    else
                        defV = 0;
                        if ismember(cn, {'TAP_MAX','TAP_MIN','K2','status'}), defV = 1; end
                        if strcmp(cn, 'SH_MIN'), defV = -360; end
                        if strcmp(cn, 'SH_MAX'), defV = 360; end
                        if ismember(cn, RIC.Properties.VariableNames)
                            TIC.(cn) = obj.evalColumn(RIC.(cn), defV);
                        elseif strcmp(cn, 'Column37') && width(RIC) >= 38
                            TIC.(cn) = obj.evalColumn(RIC{:, 38}, 0);
                        else
                            TIC.(cn) = repmat(defV, nIC, 1);
                        end
                    end
                end
                obj.tables.IntConv = TIC;
            end

            % 6. Apply .dyr model/type/MBASE overrides to Gen and IntConv
            if ~isempty(obj.dyrData) && isstruct(obj.dyrData) && isfield(obj.dyrData, 'records')
                recs = obj.dyrData.records;
                if ~isempty(obj.tables.Gen)
                    for r = 1:height(obj.tables.Gen)
                        bId = round(obj.tables.Gen.bus(r));
                        drec = VFlexPSession.findBusRecord(recs, bId, {'machine', 'device'});
                        if ~isempty(drec)
                            switch drec.model
                                case {'GENROU', 'GENSAL'}
                                    if ~strcmp(obj.tables.Gen.type(r), "SG")
                                        obj.tables.Gen.type(r) = "MGES";
                                    end
                                case 'VSM'
                                    obj.tables.Gen.type(r) = "VSM";
                                case 'GFL'
                                    obj.tables.Gen.type(r) = "GFL";
                                case 'DC_PES'
                                    obj.tables.Gen.type(r) = "DC_PES";
                            end
                            if isfield(drec.params, 'MBASE') && isnumeric(drec.params.MBASE) && ~isempty(drec.params.MBASE)
                                obj.tables.Gen.mBase(r) = double(drec.params.MBASE(1));
                            end
                        end
                    end
                end
                if ~isempty(obj.tables.IntConv)
                    for r = 1:height(obj.tables.IntConv)
                        fb = round(obj.tables.IntConv.fbus(r));
                        tb = round(obj.tables.IntConv.tbus(r));
                        irec = VFlexPSession.findBranchRecord(recs, fb, tb, 'device');
                        if ~isempty(irec)
                            switch irec.model
                                case 'INTCONV_GFL'
                                    obj.tables.IntConv.type(r) = "intconv_GFL";
                                case 'INTCONV_VSM'
                                    obj.tables.IntConv.type(r) = "intconv_VSM";
                                case 'INTCONV_GFM'
                                    obj.tables.IntConv.type(r) = "intconv_GFM";
                            end
                        end
                    end
                end
            end
        end

        function [TbusAC, TgenAC, TlineAC, TbusDC, TgenDC, TlineDC, Tintconv, TFb, Sb] = splitDomains(obj)
            %SPLITDOMAINS Separate AC and DC tables cleanly from evaluated session.tables.
            if isempty(obj.tables) || ~isfield(obj.tables, 'Bus')
                obj.evalTables();
            end
            Sb = obj.tables.Sb;
            TFb = obj.tables.Fb;
            TBusAll = obj.tables.Bus;
            isDCBus = startsWith(strtrim(TBusAll.Complexity_bus), "DC", 'IgnoreCase', true);
            TbusDC = TBusAll(isDCBus, :);
            TbusAC = TBusAll(~isDCBus, :);

            TLineAll = obj.tables.Line;
            if isempty(TLineAll) || height(TLineAll) == 0
                TlineAC = TLineAll;
                TlineDC = TLineAll;
            else
                isDCLine = startsWith(strtrim(TLineAll.Complexity), "DC", 'IgnoreCase', true);
                TlineDC = TLineAll(isDCLine, :);
                TlineAC = TLineAll(~isDCLine, :);
            end

            TGenAll = obj.tables.Gen;
            if isempty(TGenAll) || height(TGenAll) == 0
                TgenAC = TGenAll;
                TgenDC = TGenAll;
            else
                isDCGen = startsWith(strtrim(TGenAll.type), "DC", 'IgnoreCase', true);
                TgenDC = TGenAll(isDCGen, :);
                TgenAC = TGenAll(~isDCGen, :);
            end

            Tintconv = obj.tables.IntConv;
        end

        function addEvent(obj, tEv, targetType, targetId, channel, magnitude, duration, enabled)
            %ADDEVENT Append a temporal disturbance event row.
            if nargin < 7 || isempty(duration), duration = Inf; end
            if nargin < 8 || isempty(enabled), enabled = true; end
            newRow = table(logical(enabled), double(tEv), string(targetType), string(targetId), ...
                string(channel), double(magnitude), double(duration), ...
                'VariableNames', {'Enabled', 'Time', 'TargetType', 'TargetId', 'Channel', 'Magnitude', 'Duration'});
            obj.eventList = [obj.eventList; newRow];
        end

        function clearEvents(obj)
            %CLEAREVENTS Remove all temporal disturbance events.
            obj.eventList = VFlexPSession.emptyEventsTable();
        end

        function addEventExpr(obj, exprStr, duration, enabled)
            %ADDEVENTEXPR Parse and append a VFlexPy disturbance expression "t=<time> : <target> <op> <value>".
            if nargin < 3 || isempty(duration), duration = Inf; end
            if nargin < 4 || isempty(enabled), enabled = true; end
            txt = strtrim(char(string(exprStr)));
            if isempty(txt), return; end
            pat = ['^\s*t\s*=\s*([-+0-9.eE]+)\s*:\s*(\S+)\s*' ...
                   '(\*=|\+=|-=|(?<![*+\-])=)\s*' ...
                   '([-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\s*$'];
            tok = regexp(txt, pat, 'tokens', 'once');
            if isempty(tok)
                error('VFlexP:Session:InvalidEventExpr', ...
                    'Cannot parse event expression "%s". Expected "t=<time> : <target> <op> <value>".', txt);
            end
            tEv = str2double(tok{1});
            tgt = tok{2};
            op  = tok{3};
            val = str2double(tok{4});
            obj.addEvent(tEv, 'vflexpy', tgt, op, val, duration, enabled);
        end

        function setEventsFromExprs(obj, exprList)
            %SETEVENTSFROMEXPRS Replace temporal events with a list of VFlexPy event expressions.
            obj.clearEvents();
            if isempty(exprList), return; end
            strs = reshape(string(exprList), [], 1);
            for k = 1:numel(strs)
                if ~ismissing(strs(k)) && strlength(strtrim(strs(k))) > 0
                    obj.addEventExpr(strs(k));
                end
            end
        end

        function clone = cloneForSweep(obj)
            %CLONEFORSWEEP Create an independent copy of the session configuration & tables.
            clone = VFlexPSession();
            clone.casePath = obj.casePath;
            clone.caseName = obj.caseName;
            clone.rawTables = obj.rawTables;
            clone.tables = obj.tables;
            clone.variables = obj.variables;
            clone.variableRanges = obj.variableRanges;
            clone.dyrData = obj.dyrData;
            clone.unitOverrides = obj.unitOverrides;
            clone.config = obj.config;
            clone.eventList = obj.eventList;
            clone.pfResult = obj.pfResult;
        end

        function val = evalScalar(obj, rawVal, defaultVal)
            %EVALSCALAR Evaluate a numeric, string, or $variable scalar token.
            if iscell(rawVal)
                if isempty(rawVal), val = defaultVal; return; end
                rawVal = rawVal{1};
            end
            if isnumeric(rawVal)
                if isempty(rawVal) || isnan(rawVal(1))
                    val = defaultVal;
                else
                    val = double(rawVal(1));
                end
                return;
            end
            strVal = string(rawVal);
            if isempty(strVal) || ismissing(strVal(1))
                val = defaultVal;
                return;
            end
            s = strtrim(char(strVal(1)));
            if isempty(s) || strcmpi(s, '<missing>')
                val = defaultVal;
                return;
            end
            sClean = regexprep(s, '^\$', '');
            num = str2double(sClean);
            if ~isnan(num)
                val = num;
                return;
            end
            if isfield(obj.variables, sClean)
                val = double(obj.variables.(sClean));
                return;
            end
            % Substitute any $var occurrences inside an expression like "2*$H_vsm"
            expr = regexprep(s, '\$([A-Za-z_][A-Za-z0-9_]*)', '$1');
            try
                val = double(evalin('base', expr));
                if isempty(val) || ~isscalar(val) || isnan(val)
                    val = defaultVal;
                end
            catch
                val = defaultVal;
            end
        end

        function col = evalColumn(obj, rawCol, defaultVal)
            %EVALCOLUMN Evaluate a table column into a double column vector.
            if isnumeric(rawCol)
                col = double(rawCol(:));
                col(isnan(col)) = defaultVal;
                return;
            end
            n = numel(rawCol);
            col = zeros(n, 1);
            for k = 1:n
                if iscell(rawCol)
                    col(k) = obj.evalScalar(rawCol{k}, defaultVal);
                else
                    col(k) = obj.evalScalar(rawCol(k), defaultVal);
                end
            end
        end
    end

    methods (Static)
        function cols = busColumnNames()
            cols = {'bus_i','type','Pd','Qd','Gs','Bs','area','Vm','Va','baseKV','zone','Vmax','Vmin', ...
                    'Complexity_bus','Complexity_load','Latitude','Longitude'};
        end

        function cols = lineColumnNames()
            cols = {'fbus','tbus','r','x','b','rateA','rateB','rateC','ratio','angle','status','angmin','angmax','Complexity'};
        end

        function cols = genColumnNames()
            cols = {'bus','type','Pg','Qg','Qmax','Qmin','Vg','mBase','status','Pmax','Pmin', ...
                    'Pc1','Pc2','Qc1min','Qc1max','Qc2min','Qc2max','ramp_agc','ramp_10','ramp_30','ramp_q','apf'};
        end

        function cols = intConvColumnNames()
            cols = {'type','fbus','tbus','r','x','b','rateA','rateB','rateC','ratio','angle','status','angmin','angmax', ...
                    'PF','QF','PT','QT','MU_SF','MU_ST','MU_ANGMIN','MU_ANGMAX','VF_SET','VT_SET','TAP_MAX','TAP_MIN', ...
                    'CONV','BEQ','K2','BEQ_MIN','BEQ_MAX','SH_MIN','SH_MAX','GSW','ALPHA1','ALPHA2','ALPHA3','Column37'};
        end

        function rt = emptyRawTables()
            t = VFlexPSession.emptyTables();
            rt = struct( ...
                'Sb', 100, ...
                'Fb', table(1, 50, 'VariableNames', {'Area', 'Fb'}), ...
                'Bus', t.Bus, ...
                'Line', t.Line, ...
                'Gen', t.Gen, ...
                'IntConv', t.IntConv, ...
                'SecCtrl', table(), ...
                'TerCtrl', table());
        end

        function t = emptyTables()
            bCols = VFlexPSession.busColumnNames();
            TBus = array2table(zeros(0, 13), 'VariableNames', bCols(1:13));
            TBus.Complexity_bus = strings(0, 1);
            TBus.Complexity_load = strings(0, 1);
            TBus.Latitude = zeros(0, 1);
            TBus.Longitude = zeros(0, 1);

            lCols = VFlexPSession.lineColumnNames();
            TLine = array2table(zeros(0, 13), 'VariableNames', lCols(1:13));
            TLine.Complexity = strings(0, 1);

            gCols = VFlexPSession.genColumnNames();
            TGen = table(zeros(0,1), strings(0,1), 'VariableNames', {'bus', 'type'});
            for k = 3:numel(gCols)
                TGen.(gCols{k}) = zeros(0, 1);
            end

            icCols = VFlexPSession.intConvColumnNames();
            TIC = table(strings(0,1), 'VariableNames', {'type'});
            for k = 2:numel(icCols)
                TIC.(icCols{k}) = zeros(0, 1);
            end

            t = struct( ...
                'Sb', 100, ...
                'Fb', table(1, 50, 'VariableNames', {'Area', 'Fb'}), ...
                'Bus', TBus, ...
                'Line', TLine, ...
                'Gen', TGen, ...
                'IntConv', TIC);
        end

        function tev = emptyEventsTable()
            tev = table(false(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
                strings(0,1), zeros(0,1), zeros(0,1), ...
                'VariableNames', {'Enabled', 'Time', 'TargetType', 'TargetId', 'Channel', 'Magnitude', 'Duration'});
        end

        function tp = normalizeGenTypes(tp)
            tp = strtrim(string(tp));
            for k = 1:numel(tp)
                u = upper(tp(k));
                if strcmp(u, "GFR")
                    tp(k) = "GFM";
                elseif strcmp(u, "GFL")
                    tp(k) = "GFL";
                elseif ismember(u, ["SM", "GENROU", "GENSAL"])
                    tp(k) = "MGES";
                end
            end
        end

        function tp = normalizeIntConvTypes(tp)
            tp = strtrim(string(tp));
            for k = 1:numel(tp)
                u = upper(tp(k));
                if ismember(u, ["DC_GFL", "INTCONV_GFL"])
                    tp(k) = "intconv_GFL";
                elseif ismember(u, ["DC_VSM", "INTCONV_VSM"])
                    tp(k) = "intconv_VSM";
                elseif ismember(u, ["DC_GFM", "INTCONV_GFM"])
                    tp(k) = "intconv_GFM";
                end
            end
        end

        function rec = findBusRecord(records, busId, roles)
            rec = [];
            keyStr = sprintf('%d', round(busId));
            for k = 1:numel(records)
                if strcmp(records(k).key, keyStr) && ismember(records(k).role, roles)
                    rec = records(k);
                    return;
                end
            end
        end

        function rec = findBranchRecord(records, fbus, tbus, roleName)
            rec = [];
            k1 = sprintf('%d-%d', round(fbus), round(tbus));
            k2 = sprintf('%d-%d', round(tbus), round(fbus));
            for k = 1:numel(records)
                if (strcmp(records(k).key, k1) || strcmp(records(k).key, k2)) && strcmp(records(k).role, roleName)
                    rec = records(k);
                    return;
                end
            end
        end
    end
end
