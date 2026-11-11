function savePath = vflexp_save_case(session, savePath)
%VFLEXP_SAVE_CASE Save session tables, $variables, and settings to an .xlsx workbook.
%
%   savePath = vflexp_save_case(session, savePath)

if nargin < 2 || isempty(savePath)
    if ~isempty(session.casePath)
        [folder, stem, ~] = fileparts(session.casePath);
        savePath = fullfile(folder, [stem '.xlsx']);
    else
        error('VFlexP:SaveCase:NoPath', 'Specify a destination .xlsx path.');
    end
end

[folder, stem, ext] = fileparts(char(string(savePath)));
if isempty(folder), folder = pwd; end
if isempty(ext), ext = '.xlsx'; end
savePath = fullfile(folder, [stem ext]);

if ~exist(folder, 'dir')
    mkdir(folder);
end

% Write Bus data, Line data, Gen data, IntConv data, System parameters, Parameter variation
TBus = local_prepare_for_excel(session.rawTables.Bus, session.tables.Bus);
writetable(TBus, savePath, 'Sheet', 'Bus data', 'WriteMode', 'overwritesheet');

if ~isempty(session.tables.Line) && height(session.tables.Line) > 0
    TLine = local_prepare_for_excel(session.rawTables.Line, session.tables.Line);
    writetable(TLine, savePath, 'Sheet', 'Line data', 'WriteMode', 'overwritesheet');
end

TGen = local_prepare_for_excel(session.rawTables.Gen, session.tables.Gen);
writetable(TGen, savePath, 'Sheet', 'Gen data', 'WriteMode', 'overwritesheet');

if ~isempty(session.tables.IntConv) && height(session.tables.IntConv) > 0
    TIC = local_prepare_for_excel(session.rawTables.IntConv, session.tables.IntConv);
    writetable(TIC, savePath, 'Sheet', 'IntConv data', 'WriteMode', 'overwritesheet');
end

TFb = session.tables.Fb;
if isempty(TFb) || height(TFb) == 0
    TFb = table(1, 50, 'VariableNames', {'Area', 'Fb'});
end
TSys = TFb;
TSys.Sb = repmat(session.tables.Sb, height(TSys), 1);
TSys = movevars(TSys, 'Sb', 'Before', 'Area');
writetable(TSys, savePath, 'Sheet', 'System parameters', 'WriteMode', 'overwritesheet');

% Ensure all session.variables entries are mirrored in session.variableRanges
varNames = fieldnames(session.variables);
for vIdx = 1:numel(varNames)
    vn = varNames{vIdx};
    vv = double(session.variables.(vn));
    if isempty(session.variableRanges) || height(session.variableRanges) == 0
        fr = 0.8 * vv; to = 1.2 * vv; if vv == 0, fr = 0; to = 1; end
        session.variableRanges = table(string(vn), vv, fr, to, 5, ...
            'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'});
    else
        pcol = regexprep(strtrim(string(session.variableRanges.Parameter)), '^\$', '');
        hit = find(strcmp(pcol, vn), 1);
        if isempty(hit)
            fr = 0.8 * vv; to = 1.2 * vv; if vv == 0, fr = 0; to = 1; end
            session.variableRanges = [session.variableRanges; table(string(vn), vv, fr, to, 5, ...
                'VariableNames', {'Parameter', 'Value', 'From', 'To', 'N_steps'})]; %#ok<AGROW>
        else
            session.variableRanges.Value(hit) = vv;
        end
    end
end
if ~isempty(session.variableRanges) && height(session.variableRanges) > 0
    writetable(session.variableRanges, savePath, 'Sheet', 'Parameter variation', 'WriteMode', 'overwritesheet');
end

session.casePath = savePath;
session.caseName = stem;
session.log('Saved case workbook to "%s".', savePath);
end

function TOut = local_prepare_for_excel(TRaw, TEval)
if isempty(TRaw) || height(TRaw) ~= height(TEval)
    TOut = TEval;
    return;
end
TOut = TEval;
rawVars = TRaw.Properties.VariableNames;
evalVars = TEval.Properties.VariableNames;
for c = 1:numel(evalVars)
    vn = evalVars{c};
    if ismember(vn, rawVars) && (isstring(TRaw.(vn)) || iscellstr(TRaw.(vn)) || iscell(TRaw.(vn)))
        sCol = strtrim(string(TRaw.(vn)));
        % Preserve symbolic $var references in raw columns
        hasSym = startsWith(sCol, "$");
        if any(hasSym) && isnumeric(TOut.(vn))
            strOut = string(TOut.(vn));
            strOut(hasSym) = sCol(hasSym);
            TOut.(vn) = strOut;
        end
    end
end
end
