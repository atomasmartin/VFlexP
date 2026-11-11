function modalRes = vflexp_modal(session, forceRelin, exportPath)
%VFLEXP_MODAL Eigenvalue, damping, frequency, and participation factor analysis.
%
%   modalRes = vflexp_modal(session, forceRelin, exportPath)
%
%   Computes right/left eigenvectors, normalized participation matrix P
%   ([nStates x nModes], columns sum to 1), mode frequencies [Hz], damping
%   ratios [%], dominant participating states, and sortable modesTable.
%   Optionally exports results to .xlsx or .mat when exportPath is given.

if nargin < 2 || isempty(forceRelin)
    forceRelin = false;
end
if nargin < 3
    exportPath = '';
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

if nStates == 0
    modesTable = table(zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        strings(0,1), zeros(0,1), ...
        'VariableNames', {'Mode','Real_1_s','Imag_rad_s','Freq_Hz','Damping_pct','DominantState','DominantPF'});
    modalRes = struct( ...
        'eigenvalues', zeros(0,1), ...
        'realPart', zeros(0,1), ...
        'imagPart', zeros(0,1), ...
        'freqHz', zeros(0,1), ...
        'dampingPct', zeros(0,1), ...
        'V', zeros(0,0), ...
        'W', zeros(0,0), ...
        'participation', zeros(0,0), ...
        'stateLabels', stateLabels, ...
        'stateMeta', stateMeta, ...
        'modesTable', modesTable, ...
        'sortOrder', zeros(0,1), ...
        'elapsed', toc(tStart));
    session.modalResult = modalRes;
    return;
end

[V, D] = eig(A);
lam    = diag(D);

% Left eigenvectors W (rows of W satisfy W * V = I)
if rcond(V) < 1e-13
    W = pinv(V);
else
    W = inv(V);
end

% Normalized participation matrix P(k, i) = |W(i,k)| * |V(k,i)| / sum_k(...)
P = abs(W.') .* abs(V);
colSums = sum(P, 1);
colSums(colSums < eps) = 1.0;
P = P ./ colSums;

sigma  = real(lam);
omega  = imag(lam);
wn     = abs(lam);
freqHz = abs(omega) / (2 * pi);

zetaPct = zeros(nStates, 1);
nzMask  = wn > 1e-12;
zetaPct(nzMask)  = -100.0 * (sigma(nzMask) ./ wn(nzMask));
zetaPct(~nzMask) = 100.0;

domState = strings(nStates, 1);
domPF    = zeros(nStates, 1);
for m = 1:nStates
    [mxVal, mxIdx] = max(P(:, m));
    domState(m) = stateLabels(mxIdx);
    domPF(m)    = mxVal;
end

% Sort modes according to session.config.modal.sort_by
sortBy  = 'damping';
minFreq = 0.0;
if isfield(session.config, 'modal')
    if isfield(session.config.modal, 'sort_by') && ~isempty(session.config.modal.sort_by)
        sortBy = lower(char(string(session.config.modal.sort_by)));
    end
    if isfield(session.config.modal, 'min_freq') && ~isempty(session.config.modal.min_freq)
        minFreq = double(session.config.modal.min_freq);
    end
end

switch sortBy
    case 'frequency'
        [~, order] = sortrows([freqHz, zetaPct], [-1, 2]);
    case 'real'
        [~, order] = sort(sigma, 'descend');
    otherwise % 'damping' (least damped / most unstable first)
        [~, order] = sortrows([zetaPct, -sigma], [1, 2]);
end

if minFreq > 0
    keepMask = freqHz(order) >= minFreq;
    order = [order(keepMask); order(~keepMask)];
end

lam    = lam(order);
sigma  = sigma(order);
omega  = omega(order);
freqHz = freqHz(order);
zetaPct = zetaPct(order);
V      = V(:, order);
W      = W(order, :);
P      = P(:, order);
domState = domState(order);
domPF    = domPF(order);

modesTable = table((1:nStates)', sigma, omega, freqHz, zetaPct, domState, domPF, ...
    'VariableNames', {'Mode','Real_1_s','Imag_rad_s','Freq_Hz','Damping_pct','DominantState','DominantPF'});

elapsed = toc(tStart);
modalRes = struct( ...
    'eigenvalues', lam, ...
    'realPart', sigma, ...
    'imagPart', omega, ...
    'freqHz', freqHz, ...
    'dampingPct', zetaPct, ...
    'V', V, ...
    'W', W, ...
    'participation', P, ...
    'stateLabels', stateLabels, ...
    'stateMeta', stateMeta, ...
    'modesTable', modesTable, ...
    'sortOrder', order, ...
    'elapsed', elapsed);

session.modalResult = modalRes;
minDamp = min(zetaPct);
maxReal = max(sigma);
session.log('Modal analysis completed in %.3f s (%d modes, max Re=%.4g 1/s, min damping=%.2f%%).', ...
    elapsed, nStates, maxReal, minDamp);

if ~isempty(exportPath)
    local_export_modal(modalRes, char(string(exportPath)), session);
end
end

function local_export_modal(modalRes, outPath, session)
[~, ~, ext] = fileparts(outPath);
if strcmpi(ext, '.mat')
    save(outPath, 'modalRes');
else
    writetable(modalRes.modesTable, outPath, 'Sheet', 'Modes');
    nM = size(modalRes.participation, 2);
    modeColNames = arrayfun(@(i) sprintf('Mode_%d', i), 1:nM, 'UniformOutput', false);
    Tpart = [table(modalRes.stateLabels, 'VariableNames', {'State'}), ...
             array2table(modalRes.participation, 'VariableNames', modeColNames)];
    writetable(Tpart, outPath, 'Sheet', 'Participation');
end
session.log('Exported modal analysis to "%s".', outPath);
end
