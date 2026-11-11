function relRes = vflexp_relevance(session, forceRelin)
%VFLEXP_RELEVANCE Balanced-realization state relevance ranking & simplification preview.
%
%   relRes = vflexp_relevance(session, forceRelin)
%
%   Implements the state relevance coefficient of Tomas-Martin et al. (2023):
%   combines input-output Hankel singular values and balanced participation
%   factors from controllability/observability Gramians with physical modal
%   participation factors to rank every state in [0, 1] (normalized to 1).

if nargin < 2 || isempty(forceRelin)
    forceRelin = false;
end

tStart = tic;
if forceRelin || isempty(session.linResult)
    vflexp_linearise(session, forceRelin);
end

linRes      = session.linResult;
A           = linRes.A;
stateLabels = linRes.stateLabels;
stateMeta   = linRes.stateMeta;
nStates     = size(A, 1);

stabMargin = 1e-3;
threshold  = 1e-3;
inChans    = {'bus_P', 'bus_Q', 'gen_P', 'gen_V'};
outChans   = {'vac', 'vdc', 'elec'};

if isfield(session.config, 'relevance')
    rc = session.config.relevance;
    if isfield(rc, 'stability_margin') && ~isempty(rc.stability_margin)
        stabMargin = max(1e-6, double(rc.stability_margin));
    end
    if isfield(rc, 'threshold') && ~isempty(rc.threshold)
        threshold = max(0, double(rc.threshold));
    end
    if isfield(rc, 'input_channels') && ~isempty(rc.input_channels)
        inChans = cellstr(string(rc.input_channels));
    end
    if isfield(rc, 'output_channels') && ~isempty(rc.output_channels)
        outChans = cellstr(string(rc.output_channels));
    end
end

B = local_select_inputs(linRes.B, linRes.inMap, inChans);
C = local_select_outputs(linRes.C, linRes.outMap, outChans);

if nStates == 0
    rankTable = table(zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
        strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
        'VariableNames', {'Rank','StateIndex','State','Domain','ElementKind', ...
                          'SubGroup','Relevance','Accumulated','BelowThreshold'});
    relRes = struct( ...
        'stateRelevance', zeros(0,1), ...
        'modeRelevance', zeros(0,1), ...
        'eigenvalues', zeros(0,1), ...
        'hankelSV', zeros(0,1), ...
        'ranking', zeros(0,1), ...
        'accumulated', zeros(0,1), ...
        'shift', 0.0, ...
        'maxReal', -Inf, ...
        'threshold', threshold, ...
        'stateLabels', stateLabels, ...
        'stateMeta', stateMeta, ...
        'rankTable', rankTable, ...
        'planTable', local_empty_plan_table(), ...
        'elapsed', toc(tStart));
    session.relevanceResult = relRes;
    return;
end

% 1. Eigendecomposition & normalized modal participation factors P [nStates x nModes]
[V, D] = eig(A);
lam    = diag(D);
if rcond(V) < 1e-13
    W = pinv(V);
else
    W = inv(V);
end
P = abs(W.') .* abs(V);
P = P ./ max(sum(P, 1), eps);

% 2. Exponentially shifted Hurwitz matrix A_sh so Lyapunov equations have PSD solutions
maxReal = max(real(lam));
shift   = max(0.0, maxReal + stabMargin);
A_sh    = A - shift * eye(nStates);

% 3. Controllability & observability Gramians: A_sh*Wc + Wc*A_sh' + B*B' = 0
Wc = local_solve_lyap(A_sh, B * B.');
Wo = local_solve_lyap(A_sh.', C.' * C);

% 4. Symmetric PSD square-root factors Wc = Lc*Lc', Wo = Lo*Lo'
Lc = local_psd_sqrt(Wc);
Lo = local_psd_sqrt(Wo);

% 5. Balanced SVD of Lo' * Lc -> Hankel singular values & balanced modal participation
[Us, Ss, Vs] = svd(Lo.' * Lc, 'econ');
hankelSV = diag(Ss);
maxH     = max([hankelSV; eps]);

if maxH <= 1e-15 || isempty(hankelSV)
    modeRel  = ones(nStates, 1);
    stateRel = sum(P, 2);
    stateRel = stateRel / max(max(stateRel), eps);
else
    tolH   = max(maxH * 1e-13, eps);
    active = hankelSV > tolH;
    if ~any(active)
        active(1) = true;
    end
    Us_a   = Us(:, active);
    Vs_a   = Vs(:, active);
    h_a    = hankelSV(active);
    rootH  = sqrt(h_a);

    % Balanced right/left modal matrices: V_bal [r x nModes], W_bal [nModes x r]
    V_bal = (Us_a.' * (Lo.' * V)) ./ rootH;
    W_bal = ((W * Lc) * Vs_a) ./ (rootH.');

    P_bal = abs(W_bal.') .* abs(V_bal);
    pSum  = sum(P_bal, 1);
    nzCol = pSum > eps;
    P_bal(:, nzCol)  = P_bal(:, nzCol) ./ pSum(nzCol);
    P_bal(:, ~nzCol) = 0.0;

    modeRel = (h_a.' * P_bal).';
    if max(modeRel) > eps
        modeRel = modeRel / max(modeRel);
    end

    stateRel = P * modeRel;
    if max(stateRel) > eps
        stateRel = stateRel / max(stateRel);
    end
end

[sortedRel, ranking] = sort(stateRel, 'descend');
totRel = sum(sortedRel);
if totRel > eps
    accumulated = cumsum(sortedRel) / totRel;
else
    accumulated = ones(nStates, 1);
end

belowMask = sortedRel < threshold;
rankTable = table((1:nStates)', ranking, stateLabels(ranking), ...
    stateMeta.Domain(ranking), stateMeta.ElementKind(ranking), ...
    stateMeta.SubGroup(ranking), sortedRel, accumulated, belowMask, ...
    'VariableNames', {'Rank','StateIndex','State','Domain','ElementKind', ...
                      'SubGroup','Relevance','Accumulated','BelowThreshold'});

planTable = vflexp_apply_relevance(session, threshold, stateRel, stateMeta, true);

elapsed = toc(tStart);
relRes = struct( ...
    'stateRelevance', stateRel, ...
    'modeRelevance', modeRel, ...
    'eigenvalues', lam, ...
    'hankelSV', hankelSV, ...
    'ranking', ranking, ...
    'accumulated', accumulated, ...
    'shift', shift, ...
    'maxReal', maxReal, ...
    'threshold', threshold, ...
    'stateLabels', stateLabels, ...
    'stateMeta', stateMeta, ...
    'rankTable', rankTable, ...
    'planTable', planTable, ...
    'nInputs', size(B, 2), ...
    'nOutputs', size(C, 1), ...
    'elapsed', elapsed);

session.relevanceResult = relRes;
nBelow = sum(stateRel < threshold);
session.log('State relevance computed in %.3f s (%d/%d states below threshold %.3g; %d simplifications available).', ...
    elapsed, nBelow, nStates, threshold, height(planTable));
end

% =========================================================================
% Local Helpers
% =========================================================================
function Bsel = local_select_inputs(B, inMap, inChans)
cols = [];
for k = 1:numel(inChans)
    ch = strtrim(inChans{k});
    if isfield(inMap, ch)
        cols = [cols, inMap.(ch)]; %#ok<AGROW>
    end
end
cols = unique(cols, 'stable');
if isempty(cols)
    Bsel = B;
else
    Bsel = B(:, cols);
end
end

function Csel = local_select_outputs(C, outMap, outChans)
rows = [];
for k = 1:numel(outChans)
    ch = strtrim(outChans{k});
    if isfield(outMap, ch)
        rows = [rows, outMap.(ch)]; %#ok<AGROW>
    end
end
rows = unique(rows, 'stable');
if isempty(rows)
    Csel = C;
else
    Csel = C(rows, :);
end
end

function W = local_solve_lyap(A, Q)
Q = 0.5 * (Q + Q.');
if norm(Q, 'fro') == 0
    W = zeros(size(A));
    return;
end
if exist('lyap', 'file') == 2
    W = lyap(A, Q);
else
    n = size(A, 1);
    K = kron(eye(n), A) + kron(A, eye(n));
    wVec = -K \ Q(:);
    W = reshape(wVec, n, n);
end
W = 0.5 * (real(W) + real(W).');
end

function L = local_psd_sqrt(W)
W = 0.5 * (real(W) + real(W).');
[U, D] = eig(W);
d = max(real(diag(D)), 0.0);
L = real(U) .* (sqrt(d).');
end

function T = local_empty_plan_table()
T = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', {'Kind','Element','OldComplexity','NewComplexity','MaxRelevance','States'});
end
