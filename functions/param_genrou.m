function circuit = param_genrou(Tpdo, Tppdo, Tpqo, Tppqo, H, D, Xd, Xq, Xpd, Xpq, Xppd, Xl, varargin)
% PARAM_GENROU Convert standard GENROU/GENSAL synchronous machine data into
% coupled-circuit per-unit parameters for VFlexP (MGES / SG).
%
% Supports both:
%   - 'classical' conversion (textbook one-circuit-at-a-time approximation,
%     PSS/E GENROU tradition; default)
%   - 'exact' conversion (IEEE Std 1110 / Kundur 4.2 two-circuit operational
%     inductance fit matching open- and short-circuit time constants and X'')
% as well as explicit rotor circuit overrides (Lfd, Rfd, Lkd1, Rkd1, Lkq1,
% Rkq1, Lkq2, Rkq2).

p = inputParser;
addOptional(p, 'Xppq', [], @(x) isnumeric(x));
addOptional(p, 'Ra', 0.0, @(x) isnumeric(x));
addOptional(p, 'psi_mag', 0.0, @(x) isnumeric(x));
addOptional(p, 'f_base', 50.0, @(x) isnumeric(x));
addOptional(p, 'conversion', 'classical', @(x) ischar(x) || isstring(x));
addOptional(p, 'overrides', struct(), @(x) isstruct(x));
parse(p, varargin{:});

opts = p.Results;
if isempty(opts.Xppq)
    Xppq = Xppd;
else
    Xppq = opts.Xppq;
end

n = max([numel(Tpdo), numel(Tppdo), numel(Tpqo), numel(Tppqo), ...
    numel(H), numel(D), numel(Xd), numel(Xq), numel(Xpd), numel(Xpq), ...
    numel(Xppd), numel(Xl)]);

Tpdo    = local_row(Tpdo, n);
Tppdo   = local_row(Tppdo, n);
Tpqo    = local_row(Tpqo, n);
Tppqo   = local_row(Tppqo, n);
H       = local_row(H, n);
D       = local_row(D, n);
Xd      = local_row(Xd, n);
Xq      = local_row(Xq, n);
Xpd     = local_row(Xpd, n);
Xpq     = local_row(Xpq, n);
Xppd    = local_row(Xppd, n);
Xl      = local_row(Xl, n);
Xppq    = local_row(Xppq, n);
Ra      = local_row(opts.Ra, n);
psi_mag = local_row(opts.psi_mag, n);
f_base  = local_row(opts.f_base, n);

wb = 2.0 * pi .* f_base;

Lmd = Xd - Xl;
Lmq = Xq - Xl;
Lssd = Lmd + Xl;
Lssq = Lmq + Xl;

Lfd  = ((Xpd - Xl) .* (Xd - Xl)) ./ (Xd - Xpd);
Lkd1 = ((Xppd - Xl) .* (Xpd - Xl)) ./ (Xpd - Xppd);
Lkq1 = ((Xpq - Xl) .* (Xq - Xl)) ./ (Xq - Xpq);
Lkq2 = ((Xppq - Xl) .* (Xpq - Xl)) ./ (Xpq - Xppq);

conv_mode = lower(strtrim(char(opts.conversion)));
switch conv_mode
    case 'classical'
        Rfd  = (Lmd + Lfd) ./ (Tpdo .* wb);
        Rkd1 = (Lkd1 + (Lmd .* Lfd) ./ (Lmd + Lfd)) ./ (Tppdo .* wb);
        Rkq1 = (Lkq1 + Lmq) ./ (Tpqo .* wb);
        Rkq2 = (Lkq2 + (Lmq .* Lkq1) ./ (Lmq + Lkq1)) ./ (Tppqo .* wb);
    case 'exact'
        [Lfd, Lkd1, Rfd, Rkd1]   = local_exact_rotor_circuit(Xd, Xpd, Xppd, Xl, Tpdo, Tppdo, wb, 'd');
        [Lkq1, Lkq2, Rkq1, Rkq2] = local_exact_rotor_circuit(Xq, Xpq, Xppq, Xl, Tpqo, Tppqo, wb, 'q');
    otherwise
        error('VFlexP:param_genrou:InvalidConversion', ...
            'Unknown parameter conversion "%s". Choose "classical" or "exact".', conv_mode);
end

% Apply explicit equivalent-circuit overrides if provided
if isstruct(opts.overrides) && ~isempty(fieldnames(opts.overrides))
    ov = opts.overrides;
    if isfield(ov, 'Lfd')  && ~isempty(ov.Lfd),  Lfd  = local_row(ov.Lfd, n);  end
    if isfield(ov, 'Rfd')  && ~isempty(ov.Rfd),  Rfd  = local_row(ov.Rfd, n);  end
    if isfield(ov, 'Lkd1') && ~isempty(ov.Lkd1), Lkd1 = local_row(ov.Lkd1, n); end
    if isfield(ov, 'Rkd1') && ~isempty(ov.Rkd1), Rkd1 = local_row(ov.Rkd1, n); end
    if isfield(ov, 'Lkq1') && ~isempty(ov.Lkq1), Lkq1 = local_row(ov.Lkq1, n); end
    if isfield(ov, 'Rkq1') && ~isempty(ov.Rkq1), Rkq1 = local_row(ov.Rkq1, n); end
    if isfield(ov, 'Lkq2') && ~isempty(ov.Lkq2), Lkq2 = local_row(ov.Lkq2, n); end
    if isfield(ov, 'Rkq2') && ~isempty(ov.Rkq2), Rkq2 = local_row(ov.Rkq2, n); end
end

OPEN_L = 1.0e6;
OPEN_R = 1.0e5;
Lkd2 = OPEN_L * ones(1, n);
Rkd2 = OPEN_R * ones(1, n);
Lkq3 = OPEN_L * ones(1, n);
Rkq3 = OPEN_R * ones(1, n);

circuit = struct();
circuit.n        = n;
circuit.FBASE    = f_base;
circuit.H        = H;
circuit.D        = D;
circuit.Psimag_pu = psi_mag;
circuit.Lssd_pu  = Lssd;
circuit.Lssq_pu  = Lssq;
circuit.Rs_pu    = Ra;
circuit.Lmd_pu   = Lmd;
circuit.Lmq_pu   = Lmq;
circuit.Rfd_pu   = Rfd;
circuit.Lfd_pu   = Lfd;
circuit.Rkd1_pu  = Rkd1;
circuit.Lkd1_pu  = Lkd1;
circuit.Rkd2_pu  = Rkd2;
circuit.Lkd2_pu  = Lkd2;
circuit.Rkq1_pu  = Rkq1;
circuit.Lkq1_pu  = Lkq1;
circuit.Rkq2_pu  = Rkq2;
circuit.Lkq2_pu  = Lkq2;
circuit.Rkq3_pu  = Rkq3;
circuit.Lkq3_pu  = Lkq3;
circuit.avr_to_vfd = Rfd ./ Lmd;
circuit.conversion = upper(conv_mode);
end

function v = local_row(x, n)
x = reshape(double(x), 1, []);
if isscalar(x)
    v = repmat(x, 1, n);
elseif numel(x) == n
    v = x;
else
    error('VFlexP:param_genrou:DimensionMismatch', ...
        'Expected scalar or vector of length %d, got length %d.', n, numel(x));
end
end

function [L1, L2, r1_pu, r2_pu] = local_exact_rotor_circuit(X, Xp, Xpp, Xl, T1o, T2o, wb, axis_name)
% Exact two-circuit rotor ladder matching operational inductance poles
% (T'o, T''o) and short-circuit zeros (T' = T'o*X'/X, T'' = T''o*X''/X').
La = X - Xl;
T1 = T1o .* Xp ./ X;
T2 = T2o .* Xpp ./ Xp;
a1 = (T1o + T2o) - (T1 + T2);
a2 = T1o .* T2o - T1 .* T2;
b1 = X .* (T1 + T2) - Xl .* (T1o + T2o);
b2 = X .* T1 .* T2 - Xl .* T1o .* T2o;
disc = b1 .* b1 - 4.0 .* b2 .* La;
if any(disc < 0 | a1 <= 0 | a2 <= 0 | b2 <= 0)
    error('VFlexP:param_genrou:ExactConversionFailed', ...
        '%s-axis data cannot be met exactly by a two-circuit rotor; check data or use classical conversion.', axis_name);
end
root = sqrt(disc);
tau1 = 2.0 .* b2 ./ (b1 - root);
tau2 = 2.0 .* b2 ./ (b1 + root);
q = a2 ./ a1;
rho = (q - tau2) ./ (tau1 - q);
r2 = La .* La .* (rho + 1.0) ./ (X .* a1 .* rho);
r1 = rho .* r2;
L1 = r1 .* tau1;
L2 = r2 .* tau2;
if any(L1 <= 0 | L2 <= 0 | r1 <= 0 | r2 <= 0)
    error('VFlexP:param_genrou:NegativeRotorParameter', ...
        'Exact %s-axis conversion gives negative inductance or resistance; use classical conversion.', axis_name);
end
r1_pu = r1 ./ wb;
r2_pu = r2 ./ wb;
end

