function sweepRes = vflexp_sweep(session, sweepSpec, progressCb)
%VFLEXP_SWEEP Parametric modal or state-relevance sweep (dynamic vs static PF).
%
%   sweepRes = vflexp_sweep(session, sweepSpec, progressCb)
%
%   sweepSpec fields:
%     kind       : 'modal' (default) | 'relevance'
%     mode       : 'keep'  (keep current PF operating point — dynamic sweep) |
%                  'recompute' (re-run PF at every step — static+dynamic sweep)
%     targetType : 'variable' | 'unitParam' | 'loadScale' | 'tableCell'
%     paramName  : e.g. 'H_vsm' (for variable), 'VSM.Hv' (for unitParam),
%                  'Pd_Qd' (for loadScale), or 'Gen.Pg(1)' (for tableCell)
%     values     : numeric vector of parameter values (or from/to/nSteps)
%     from, to, nSteps : used if values is omitted

if nargin < 3
    progressCb = [];
end

tStart = tic;
kind = 'modal';
if isfield(sweepSpec, 'kind') && ~isempty(sweepSpec.kind)
    kind = lower(char(string(sweepSpec.kind)));
end

pfMode = 'keep';
if isfield(sweepSpec, 'mode') && ~isempty(sweepSpec.mode)
    pfMode = lower(char(string(sweepSpec.mode)));
end

targetType = 'variable';
if isfield(sweepSpec, 'targetType') && ~isempty(sweepSpec.targetType)
    targetType = char(string(sweepSpec.targetType));
end

paramName = '';
if isfield(sweepSpec, 'paramName')
    paramName = strtrim(char(string(sweepSpec.paramName)));
end

if isfield(sweepSpec, 'values') && ~isempty(sweepSpec.values)
    vals = reshape(double(sweepSpec.values), 1, []);
else
    fromV  = double(sweepSpec.from);
    toV    = double(sweepSpec.to);
    nSteps = max(2, round(double(sweepSpec.nSteps)));
    vals   = linspace(fromV, toV, nSteps);
end
nSteps = numel(vals);

% Ensure baseline PF is available when mode == 'keep'
if strcmp(pfMode, 'keep') && (isempty(session.pfResult) || ~session.pfResult.converged)
    vflexp_run_pf(session);
end

% Work on an isolated clone so the user's active session tables/variables are untouched
restoreBaseVars = onCleanup(@() session.syncVariablesToBase());
if isempty(restoreBaseVars), return; end
workSess = session.cloneForSweep();
baseBusPd = workSess.tables.Bus.Pd;
baseBusQd = workSess.tables.Bus.Qd;

stepResults = cell(1, nSteps);
converged   = true(1, nSteps);
stepErrors  = strings(1, nSteps);

for s = 1:nSteps
    v = vals(s);
    if ~isempty(progressCb)
        progressCb(s, nSteps, v);
    end

    try
        local_apply_step_value(workSess, targetType, paramName, v, baseBusPd, baseBusQd);
        if strcmp(pfMode, 'recompute')
            pfR = vflexp_run_pf(workSess);
            if ~pfR.converged
                converged(s)  = false;
                stepErrors(s) = "Power Flow did not converge";
                continue;
            end
        end

        vflexp_build_and_init(workSess);
        vflexp_linearise(workSess, false);

        if strcmp(kind, 'relevance')
            stepResults{s} = vflexp_relevance(workSess, false);
        else
            stepResults{s} = vflexp_modal(workSess, false);
        end
    catch ME
        converged(s)  = false;
        stepErrors(s) = string(ME.message);
    end
end

% Assemble trajectory matrices across converged steps
firstIdx = find(converged & ~cellfun(@isempty, stepResults), 1);
if isempty(firstIdx)
    error('VFlexP:Sweep:AllStepsFailed', ...
        'All %d sweep steps failed. First error: %s', nSteps, stepErrors(1));
end

refRes      = stepResults{firstIdx};
stateLabels = refRes.stateLabels;
stateMeta   = refRes.stateMeta;
nStates     = numel(stateLabels);

elapsed = toc(tStart);
sweepRes = struct( ...
    'kind', kind, ...
    'mode', pfMode, ...
    'targetType', targetType, ...
    'paramName', paramName, ...
    'values', vals, ...
    'converged', converged, ...
    'errors', stepErrors, ...
    'stateLabels', stateLabels, ...
    'stateMeta', stateMeta, ...
    'stepResults', {stepResults}, ...
    'elapsed', elapsed);

if strcmp(kind, 'relevance')
    relMat   = nan(nStates, nSteps);
    shiftVec = nan(1, nSteps);
    for s = 1:nSteps
        if converged(s) && ~isempty(stepResults{s}) && numel(stepResults{s}.stateRelevance) == nStates
            relMat(:, s) = stepResults{s}.stateRelevance;
            shiftVec(s)  = stepResults{s}.shift;
        end
    end
    sweepRes.relevanceMatrix = relMat;
    sweepRes.shifts          = shiftVec;
    if strcmp(pfMode, 'recompute')
        session.sweepRelStatic = sweepRes;
    else
        session.sweepRelDynamic = sweepRes;
    end
else
    eigMat  = complex(nan(nStates, nSteps), nan(nStates, nSteps));
    dampMat = nan(nStates, nSteps);
    freqMat = nan(nStates, nSteps);
    prevEig = [];
    for s = 1:nSteps
        if converged(s) && ~isempty(stepResults{s}) && numel(stepResults{s}.eigenvalues) == nStates
            curEig = stepResults{s}.eigenvalues(:);
            if isempty(prevEig)
                ord = (1:nStates)';
            else
                ord = local_match_eigenvalues(prevEig, curEig);
            end
            matchedEig    = curEig(ord);
            eigMat(:, s)  = matchedEig;
            freqMat(:, s) = abs(imag(matchedEig)) / (2 * pi);
            wn            = abs(matchedEig);
            z             = 100.0 * ones(nStates, 1);
            nz            = wn > 1e-12;
            z(nz)         = -100.0 * real(matchedEig(nz)) ./ wn(nz);
            dampMat(:, s) = z;
            prevEig       = matchedEig;
        end
    end
    sweepRes.eigenvalueTraj = eigMat;
    sweepRes.dampingTraj    = dampMat;
    sweepRes.freqTraj       = freqMat;
    if strcmp(pfMode, 'recompute')
        session.sweepModalStatic = sweepRes;
    else
        session.sweepModalDynamic = sweepRes;
    end
end

session.log('Completed %s sweep (%s, %s, %d/%d steps converged) in %.3f s.', ...
    kind, pfMode, paramName, sum(converged), nSteps, elapsed);
end

% =========================================================================
% Local Helpers
% =========================================================================
function local_apply_step_value(sess, targetType, paramName, val, basePd, baseQd)
switch lower(targetType)
    case 'variable'
        sess.setVariable(paramName, val);

    case 'unitparam'
        % Format: 'UnitName.FieldName' or 'UnitName.FieldName(idx)'
        tok = regexp(paramName, '^([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)(?:\((\d+)\))?$', 'tokens', 'once');
        if isempty(tok)
            error('VFlexP:Sweep:InvalidUnitParam', ...
                'Expected "Unit.Param" or "Unit.Param(idx)", got "%s".', paramName);
        end
        uName = tok{1};
        fName = tok{2};
        if ~isfield(sess.unitOverrides, uName) || ~isstruct(sess.unitOverrides.(uName))
            sess.unitOverrides.(uName) = struct();
        end
        if numel(tok) < 3 || isempty(tok{3})
            sess.unitOverrides.(uName).(fName) = val;
        else
            idx = str2double(tok{3});
            if isfield(sess.unitOverrides.(uName), fName)
                vec = sess.unitOverrides.(uName).(fName);
            else
                vec = val * ones(1, idx);
            end
            vec(idx) = val;
            sess.unitOverrides.(uName).(fName) = vec;
        end
        sess.invalidateDynamicResults();

    case 'loadscale'
        sess.tables.Bus.Pd = basePd * val;
        sess.tables.Bus.Qd = baseQd * val;
        sess.rawTables.Bus.Pd = sess.tables.Bus.Pd;
        sess.rawTables.Bus.Qd = sess.tables.Bus.Qd;
        sess.invalidateDynamicResults();

    case 'tablecell'
        % Format: 'TableName.ColName(row)' e.g. 'Gen.Pg(1)' or 'Line.x(2)'
        tok = regexp(paramName, '^([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\((\d+)\)$', 'tokens', 'once');
        if isempty(tok)
            error('VFlexP:Sweep:InvalidTableCell', ...
                'Expected "Table.Col(row)", got "%s".', paramName);
        end
        tName = tok{1};
        cName = tok{2};
        rIdx  = str2double(tok{3});
        sess.tables.(tName).(cName)(rIdx) = val;
        if isfield(sess.rawTables, tName) && ismember(cName, sess.rawTables.(tName).Properties.VariableNames)
            if iscell(sess.rawTables.(tName).(cName))
                sess.rawTables.(tName).(cName){rIdx} = val;
            elseif isstring(sess.rawTables.(tName).(cName))
                sess.rawTables.(tName).(cName)(rIdx) = string(num2str(val, '%.12g'));
            else
                sess.rawTables.(tName).(cName)(rIdx) = val;
            end
        end
        sess.invalidateDynamicResults();

    otherwise
        error('VFlexP:Sweep:UnknownTargetType', 'Unknown sweep targetType: "%s".', targetType);
end
end

function ord = local_match_eigenvalues(prevEig, curEig)
% Greedy minimum-distance assignment in complex plane to keep mode loci smooth
n   = numel(prevEig);
D   = abs(prevEig(:) - (curEig(:).'));
ord = zeros(n, 1);
usedCols = false(1, n);
usedRows = false(n, 1);

% Sort all pairs by distance
[rGrid, cGrid] = ndgrid(1:n, 1:n);
[~, sIdx] = sort(D(:), 'ascend');
assigned = 0;
for k = 1:numel(sIdx)
    r = rGrid(sIdx(k));
    c = cGrid(sIdx(k));
    if ~usedRows(r) && ~usedCols(c)
        ord(r)      = c;
        usedRows(r) = true;
        usedCols(c) = true;
        assigned    = assigned + 1;
        if assigned == n
            break;
        end
    end
end
end
