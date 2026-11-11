function session = vflexp_load_case(sessionOrPath, casePath)
%VFLEXP_LOAD_CASE Load a VFlexP case (.xlsx + optional .dyr + optional .json, or MATPOWER .m).
%
%   session = vflexp_load_case(casePath)
%   session = vflexp_load_case(session, casePath)

if nargin == 1
    casePath = sessionOrPath;
    session = VFlexPSession();
else
    session = sessionOrPath;
    if isempty(session) || ~isa(session, 'VFlexPSession')
        session = VFlexPSession();
    end
end

fullPath = char(string(casePath));
[folder, stem, ext] = fileparts(fullPath);
if isempty(folder)
    if isfile(fullPath)
        [folder, stem, ext] = fileparts(which(fullPath));
        if isempty(folder)
            folder = pwd;
        end
    else
        folder = pwd;
    end
end
if isempty(ext)
    ext = '.xlsx';
end

session.resetDefaults();
session.caseName = stem;

if strcmpi(ext, '.m')
    mFile = fullfile(folder, [stem '.m']);
    if ~isfile(mFile)
        error('VFlexP:LoadCase:FileNotFound', 'MATPOWER case file not found: %s', mFile);
    end
    session.casePath = mFile;
    oldDir = pwd;
    cleanObj = onCleanup(@() cd(oldDir));
    cd(folder);
    mpc = feval(stem);
    clear cleanObj;
    local_load_from_mpc(session, mpc);
    session.evalTables();
    session.log('Loaded MATPOWER case "%s" (%d buses, %d generators).', ...
        stem, height(session.tables.Bus), height(session.tables.Gen));
    return;
end

xlsxFile = fullfile(folder, [stem '.xlsx']);
dyrFile  = fullfile(folder, [stem '.dyr']);
jsonFile = fullfile(folder, [stem '.json']);

if ~isfile(xlsxFile)
    error('VFlexP:LoadCase:MissingExcel', 'Workbook file not found: %s', xlsxFile);
end
session.casePath = xlsxFile;

% 1. Optional sibling <stem>.json configuration
if isfile(jsonFile)
    js = fun_read_case_json(jsonFile);
    if ~isempty(js.t_end)
        session.config.sim.t_end = js.t_end;
    end
    session.config.sim.method   = js.solver.method;
    session.config.sim.rtol     = js.solver.rtol;
    session.config.sim.atol     = js.solver.atol;
    session.config.sim.max_step = js.solver.max_step;
    session.config.powerflow.tolerance        = js.powerflow.tolerance;
    session.config.powerflow.max_iterations   = js.powerflow.max_iterations;
    session.config.powerflow.enforce_q_limits = js.powerflow.enforce_q_limits;
end

% 2. Optional sibling <stem>.dyr dynamic data & $definitions
if isfile(dyrFile)
    session.dyrData = fun_read_dyr(dyrFile);
    for k = 1:numel(session.dyrData.definition_names)
        vname = session.dyrData.definition_names{k};
        vval  = session.dyrData.definitions.(vname);
        session.variables.(vname) = vval;
    end
else
    session.dyrData = [];
end

% 3. Read Excel sheets preserving $variable strings
shNames = sheetnames(xlsxFile);

if ismember("Parameter variation", shNames)
    TParam = local_read_sheet_raw(xlsxFile, "Parameter variation");
    nP = height(TParam);
    pNames = strings(nP, 1);
    pVals  = zeros(nP, 1);
    pFrom  = zeros(nP, 1);
    pTo    = zeros(nP, 1);
    pSteps = repmat(5, nP, 1);
    for r = 1:nP
        nm = regexprep(strtrim(char(string(TParam{r, 1}))), '^\$', '');
        pNames(r) = string(nm);
        vl = session.evalScalar(TParam{r, 2}, 0);
        pVals(r) = vl;
        if width(TParam) >= 3, pFrom(r) = session.evalScalar(TParam{r, 3}, 0.8*vl); end
        if width(TParam) >= 4, pTo(r)   = session.evalScalar(TParam{r, 4}, 1.2*vl); end
        if width(TParam) >= 5, pSteps(r) = max(2, round(session.evalScalar(TParam{r, 5}, 5))); end
        if ~isempty(nm) && isvarname(nm)
            session.variables.(nm) = vl;
        end
    end
    validRows = pNames ~= "";
    session.variableRanges = table(pNames(validRows), pVals(validRows), pFrom(validRows), ...
        pTo(validRows), pSteps(validRows), ...
        'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
end

% Append any .dyr $definitions not already in variableRanges
if ~isempty(session.dyrData) && ~isempty(session.dyrData.definition_names)
    existing = string.empty(0, 1);
    if ~isempty(session.variableRanges)
        existing = regexprep(strtrim(string(session.variableRanges.Parameter)), '^\$', '');
    end
    for k = 1:numel(session.dyrData.definition_names)
        nm = session.dyrData.definition_names{k};
        if ~ismember(string(nm), existing)
            vl = session.dyrData.definitions.(nm);
            fr = 0.8 * vl;
            to = 1.2 * vl;
            if vl == 0, fr = 0.0; to = 1.0; end
            addRow = table(string(nm), vl, fr, to, 5, ...
                'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
            session.variableRanges = [session.variableRanges; addRow];
        end
    end
end

% Apply optional variables & sim overrides from sibling <stem>.json
if isfile(jsonFile)
    rawJson = jsondecode(fileread(jsonFile));
    if isfield(rawJson, 'sim') && isstruct(rawJson.sim)
        if isfield(rawJson.sim, 't_end') && ~isempty(rawJson.sim.t_end)
            session.config.sim.t_end = double(rawJson.sim.t_end);
        end
        if isfield(rawJson.sim, 'method') && ~isempty(rawJson.sim.method)
            session.config.sim.method = char(string(rawJson.sim.method));
        end
    end
    if isfield(rawJson, 'variables') && isstruct(rawJson.variables)
        jVars = fieldnames(rawJson.variables);
        for jv = 1:numel(jVars)
            nm = jVars{jv};
            entry = rawJson.variables.(nm);
            if isstruct(entry) && isfield(entry, 'value')
                vl = double(entry.value);
                fr = 0.8 * vl; to = 1.2 * vl; st = 5;
                if isfield(entry, 'min'),    fr = double(entry.min); end
                if isfield(entry, 'max'),    to = double(entry.max); end
                if isfield(entry, 'points'), st = max(2, round(double(entry.points))); end
            elseif isnumeric(entry) && ~isempty(entry)
                vl = double(entry(1));
                fr = 0.8 * vl; to = 1.2 * vl; st = 5;
                if vl == 0, fr = 0.0; to = 1.0; end
            else
                continue;
            end
            session.variables.(nm) = vl;
            if isempty(session.variableRanges) || height(session.variableRanges) == 0
                session.variableRanges = table(string(nm), vl, fr, to, st, ...
                    'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
            else
                pcol = regexprep(strtrim(string(session.variableRanges.Parameter)), '^\$', '');
                hit = find(strcmp(pcol, nm), 1);
                if isempty(hit)
                    session.variableRanges = [session.variableRanges; table(string(nm), vl, fr, to, st, ...
                        'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'})]; %#ok<AGROW>
                else
                    session.variableRanges.Value(hit)   = vl;
                    session.variableRanges.From(hit)    = fr;
                    session.variableRanges.To(hit)      = to;
                    session.variableRanges.N_steps(hit) = st;
                end
            end
        end
    end
end

session.syncVariablesToBase();

% Refresh .dyr $references in case Parameter variation overrode a .dyr $var
if ~isempty(session.dyrData) && isfield(session.dyrData, 'records')
    for r = 1:numel(session.dyrData.records)
        refs = session.dyrData.records(r).references;
        if isstruct(refs)
            rf = fieldnames(refs);
            for j = 1:numel(rf)
                vn = refs.(rf{j});
                if isfield(session.variables, vn)
                    session.dyrData.records(r).params.(rf{j}) = session.variables.(vn);
                end
            end
        end
    end
end

if ismember("System parameters", shNames)
    TSys = local_read_sheet_raw(xlsxFile, "System parameters");
    if ismember('Sb', TSys.Properties.VariableNames) && height(TSys) >= 1
        session.rawTables.Sb = TSys.Sb(1);
    end
    if all(ismember({'Area', 'Fb'}, TSys.Properties.VariableNames))
        session.rawTables.Fb = TSys(:, {'Area', 'Fb'});
    end
end

session.rawTables.Bus = local_read_sheet_raw(xlsxFile, "Bus data");
if ismember("Line data", shNames)
    session.rawTables.Line = local_read_sheet_raw(xlsxFile, "Line data");
end
session.rawTables.Gen = local_read_sheet_raw(xlsxFile, "Gen data");
if ismember("IntConv data", shNames)
    session.rawTables.IntConv = local_read_sheet_raw(xlsxFile, "IntConv data");
end
if ismember("Secondary Control", shNames)
    TSec = local_read_sheet_raw(xlsxFile, "Secondary Control");
    session.rawTables.SecCtrl = TSec;
    if ismember('SC_structure', TSec.Properties.VariableNames) && height(TSec) >= 1
        session.config.sim.sc_structure = char(strtrim(string(TSec.SC_structure(1))));
    end
    if ismember('Define_SC_in_MATLAB', TSec.Properties.VariableNames) && height(TSec) >= 1
        session.config.sim.define_sc_in_matlab = logical(session.evalScalar(TSec.Define_SC_in_MATLAB(1), 0));
    end
end
if ismember("Tertiary Control", shNames)
    TTer = local_read_sheet_raw(xlsxFile, "Tertiary Control");
    session.rawTables.TerCtrl = TTer;
    if ismember('TC_structure', TTer.Properties.VariableNames) && height(TTer) >= 1
        session.config.sim.tc_structure = char(strtrim(string(TTer.TC_structure(1))));
    end
end

% Evaluate numeric tables & apply .dyr model types
session.evalTables();

% Also apply .dyr converter loss parameters (ALPHA1..3, GSW) into IntConv if present
local_apply_dyr_intconv_static(session);

extras = strings(0, 1);
if isfile(dyrFile),  extras(end+1) = ".dyr"; end
if isfile(jsonFile), extras(end+1) = ".json"; end
if isempty(extras)
    extTag = ".xlsx";
else
    extTag = ".xlsx + " + strjoin(extras, " + ");
end
session.log('Loaded case "%s" (%s): %d buses, %d lines, %d gens, %d intconvs.', ...
    stem, extTag, height(session.tables.Bus), height(session.tables.Line), ...
    height(session.tables.Gen), height(session.tables.IntConv));
end

% -------------------------------------------------------------------------
function T = local_read_sheet_raw(xlsxFile, sheetName)
opts = detectImportOptions(xlsxFile, 'Sheet', sheetName, 'PreserveVariableNames', true);
% Keep text/mixed columns as string so "$var" tokens are never lost as NaN
T = readtable(xlsxFile, opts);
% Fix column names if MATLAB modified them
for c = 1:width(T)
    if iscell(T{:, c})
        T.(T.Properties.VariableNames{c}) = string(T{:, c});
    end
end
end

function local_apply_dyr_intconv_static(session)
if isempty(session.dyrData) || isempty(session.tables.IntConv)
    return;
end
recs = session.dyrData.records;
for r = 1:height(session.tables.IntConv)
    fb = round(session.tables.IntConv.fbus(r));
    tb = round(session.tables.IntConv.tbus(r));
    irec = VFlexPSession.findBranchRecord(recs, fb, tb, 'device');
    if ~isempty(irec) && isstruct(irec.params)
        for fn = {'ALPHA1', 'ALPHA2', 'ALPHA3', 'GSW'}
            col = fn{1};
            if isfield(irec.params, col) && isnumeric(irec.params.(col))
                val = double(irec.params.(col)(1));
                session.tables.IntConv.(col)(r) = val;
                if ismember(col, session.rawTables.IntConv.Properties.VariableNames)
                    session.rawTables.IntConv.(col)(r) = val;
                end
            end
        end
    end
end
end

function local_load_from_mpc(session, mpc)
session.rawTables.Sb = mpc.baseMVA;
bCols = VFlexPSession.busColumnNames();
TBus = array2table(mpc.bus(:, 1:13), 'VariableNames', bCols(1:13));
nB = height(TBus);
TBus.Complexity_bus = repmat("static omega", nB, 1);
TBus.Complexity_load = repmat("series RL static omega", nB, 1);
TBus.Latitude = zeros(nB, 1);
TBus.Longitude = zeros(nB, 1);

gCols = VFlexPSession.genColumnNames();
nG = size(mpc.gen, 1);
nGc = min(size(mpc.gen, 2), numel(gCols) - 1);
genMat = zeros(nG, numel(gCols) - 1);
genMat(:, 1:nGc) = mpc.gen(:, 1:nGc);
TGen = table(genMat(:, 1), repmat("GFM", nG, 1), 'VariableNames', {'bus', 'type'});
for k = 3:numel(gCols)
    TGen.(gCols{k}) = genMat(:, k - 1);
end

TLine = VFlexPSession.emptyTables().Line;
TIC   = VFlexPSession.emptyTables().IntConv;
if ~isempty(mpc.branch)
    lCols = VFlexPSession.lineColumnNames();
    if size(mpc.branch, 2) > 13
        idxLines = mpc.branch(:, 26) == 0;
        idxDCBr  = (mpc.branch(:, 4) == 0) & idxLines;
        dcBuses  = unique([mpc.branch(idxDCBr, 1); mpc.branch(idxDCBr, 2)]);
        isDCBus  = ismember(TBus.bus_i, dcBuses);
        TBus.Complexity_bus(isDCBus) = "DC static";
        TBus.Complexity_load(isDCBus) = "DC ConstP";

        TLine = array2table(mpc.branch(idxLines, 1:13), 'VariableNames', lCols(1:13));
        TLine.Complexity = repmat("static omega", height(TLine), 1);
        TLine.Complexity(mpc.branch(idxLines, 4) == 0) = "DC static";

        idxIC = mpc.branch(:, 26) ~= 0;
        if any(idxIC)
            icCols = VFlexPSession.intConvColumnNames();
            nIC = sum(idxIC);
            brIC = mpc.branch(idxIC, :);
            nBrCols = min(size(brIC, 2), numel(icCols) - 1);
            icMat = zeros(nIC, numel(icCols) - 1);
            icMat(:, 1:nBrCols) = brIC(:, 1:nBrCols);
            TIC = [table(repmat("intconv_GFL", nIC, 1), 'VariableNames', {'type'}), ...
                   array2table(icMat, 'VariableNames', icCols(2:end))];
        end
    else
        TLine = array2table(mpc.branch(:, 1:13), 'VariableNames', lCols(1:13));
        TLine.Complexity = repmat("static omega", height(TLine), 1);
    end
end

areas = unique(TBus.area);
if isempty(areas), areas = 1; end
session.rawTables.Fb = table(areas(:), repmat(50, numel(areas), 1), 'VariableNames', {'Area', 'Fb'});
session.rawTables.Bus = TBus;
session.rawTables.Gen = TGen;
session.rawTables.Line = TLine;
session.rawTables.IntConv = TIC;
end
