function pfRes = vflexp_run_pf(session)
%VFLEXP_RUN_PF Execute AC or hybrid AC/DC power flow via MATPOWER / MATPOWER-FUBM.
%
%   pfRes = vflexp_run_pf(session)

tStart = tic;
session.evalTables();
mpc = local_session_to_mpc(session);

mpopt = mpoption('verbose', 0, 'out.all', 0);
if isfield(session.config, 'powerflow')
    pfc = session.config.powerflow;
    if isfield(pfc, 'tolerance') && ~isempty(pfc.tolerance)
        mpopt.pf.tol = double(pfc.tolerance);
    end
    if isfield(pfc, 'max_iterations') && ~isempty(pfc.max_iterations)
        mpopt.pf.nr.max_it = round(double(pfc.max_iterations));
    end
    if isfield(pfc, 'enforce_q_limits') && ~isempty(pfc.enforce_q_limits)
        mpopt.pf.enforce_q_lims = double(logical(pfc.enforce_q_limits));
    end
end

mpcRes = runpf(mpc, mpopt);
elapsed = toc(tStart);

pfRes = local_unpack_pf_result(session, mpc, mpcRes, elapsed, 'PF');
session.pfResult = pfRes;
session.invalidateDynamicResults();

if pfRes.converged
    session.log('Power Flow converged in %.3f s (Pgen = %.2f MW, Pload = %.2f MW, Ploss = %.3f MW).', ...
        elapsed, pfRes.summary.Pgen_MW, pfRes.summary.Pload_MW, pfRes.summary.Ploss_MW);
else
    session.log('WARNING: MATPOWER Power Flow did NOT converge.');
end
end

% -------------------------------------------------------------------------
function mpc = local_session_to_mpc(session)
TBus = session.tables.Bus;
TLine = session.tables.Line;
TGen = session.tables.Gen;
TIC = session.tables.IntConv;

isDCBus = startsWith(strtrim(TBus.Complexity_bus), "DC", 'IgnoreCase', true);
TBusMP = TBus;
TBusMP.Bs(isDCBus) = 0;

mpc = struct();
mpc.version = '2';
mpc.baseMVA = double(session.tables.Sb);
mpc.bus = table2array(TBusMP(:, 1:13));
mpc.gen = table2array(TGen(:, [1, 3:width(TGen)]));
mpc.bus = local_ensure_slack_bus(mpc.bus, mpc.gen, isDCBus, session);
mpc.gencost = repmat([2, 0, 0, 3, 0, 20, 0], height(TGen), 1);

hasLines = ~isempty(TLine) && height(TLine) > 0;
hasIC    = ~isempty(TIC) && height(TIC) > 0;

if ~hasLines && ~hasIC
    mpc.branch = [];
else
    if hasLines
        TLineMP = TLine;
        isDCLine = startsWith(strtrim(TLineMP.Complexity), "DC", 'IgnoreCase', true);
        TLineMP.x(isDCLine) = 0;
        TLineMP.b(isDCLine) = 0;
        brLines = table2array(TLineMP(:, 1:13));
    else
        brLines = zeros(0, 13);
    end

    if hasIC
        brIC = table2array(TIC(:, 2:width(TIC)));
        nFubmCols = size(brIC, 2);
        if hasLines
            padLines = zeros(size(brLines, 1), nFubmCols);
            padLines(:, 1:13) = brLines;
            % TAP_MAX=25, TAP_MIN=26? Wait: in FUBM branch (37 cols):
            % cols 24=TAP_MAX, 25=TAP_MIN, 28=K2, 31=SH_MIN, 32=SH_MAX
            icNames = TIC.Properties.VariableNames(2:end);
            idxTapMax = find(strcmp(icNames, 'TAP_MAX'), 1);
            idxTapMin = find(strcmp(icNames, 'TAP_MIN'), 1);
            idxK2     = find(strcmp(icNames, 'K2'), 1);
            idxShMin  = find(strcmp(icNames, 'SH_MIN'), 1);
            idxShMax  = find(strcmp(icNames, 'SH_MAX'), 1);
            if ~isempty(idxTapMax), padLines(:, idxTapMax) = 1; end
            if ~isempty(idxTapMin), padLines(:, idxTapMin) = 1; end
            if ~isempty(idxK2),     padLines(:, idxK2)     = 1; end
            if ~isempty(idxShMin),  padLines(:, idxShMin)  = -360; end
            if ~isempty(idxShMax),  padLines(:, idxShMax)  = 360; end
            mpc.branch = [padLines; brIC];
        else
            mpc.branch = brIC;
        end
    else
        mpc.branch = brLines;
    end
end
end

function busMat = local_ensure_slack_bus(busMat, genMat, isDCBus, session)
if isempty(busMat) || isempty(genMat), return; end
onGenBuses = genMat(genMat(:, 8) > 0, 1);
acMask = ~isDCBus(:);
hasRefOrPV = any(acMask & ismember(busMat(:, 2), [2, 3]) & ismember(busMat(:, 1), onGenBuses));
if ~hasRefOrPV
    cand = find(acMask & ismember(busMat(:, 1), onGenBuses), 1);
    if ~isempty(cand)
        busMat(cand, 2) = 3;
        session.log('Auto-assigned AC Bus %d as Slack (type=3) for Power Flow.', round(busMat(cand, 1)));
    end
end
end

% -------------------------------------------------------------------------
function pfRes = local_unpack_pf_result(session, mpcIn, mpcRes, elapsed, solverKind)
Sb = mpcIn.baseMVA;
converged = isfield(mpcRes, 'success') && mpcRes.success == 1;

brComplete = mpcRes.branch;
intconv_PF = struct('include', false, 'vmod', [], 'vangle', [], 'P', [], 'Q', [], 'Pdc', [], 'from', [], 'to', []);

if converged
    % Update Bus Vm and Va (preserving DC Bs)
    session.tables.Bus.Vm = mpcRes.bus(:, 8);
    session.tables.Bus.Va = mpcRes.bus(:, 9);
    session.rawTables.Bus = local_update_raw_col(session.rawTables.Bus, 'Vm', mpcRes.bus(:, 8));
    session.rawTables.Bus = local_update_raw_col(session.rawTables.Bus, 'Va', mpcRes.bus(:, 9));

    % Update Gen Pg, Qg, Vg
    session.tables.Gen.Pg = mpcRes.gen(:, 2);
    session.tables.Gen.Qg = mpcRes.gen(:, 3);
    session.tables.Gen.Vg = mpcRes.gen(:, 6);
    session.rawTables.Gen = local_update_raw_col(session.rawTables.Gen, 'Pg', mpcRes.gen(:, 2));
    session.rawTables.Gen = local_update_raw_col(session.rawTables.Gen, 'Qg', mpcRes.gen(:, 3));
    session.rawTables.Gen = local_update_raw_col(session.rawTables.Gen, 'Vg', mpcRes.gen(:, 6));

    if ~isempty(brComplete) && size(brComplete, 2) >= 37
        idxIC = find(brComplete(:, 26) ~= 0);
        nIC = numel(idxIC);
        if nIC > 0
            intconv_PF.include = true;
            intconv_PF.vmod   = zeros(1, nIC);
            intconv_PF.vangle = zeros(1, nIC);
            intconv_PF.P      = zeros(1, nIC);
            intconv_PF.Q      = zeros(1, nIC);
            intconv_PF.Pdc    = zeros(1, nIC);
            intconv_PF.from   = brComplete(idxIC, 1)';
            intconv_PF.to     = brComplete(idxIC, 2)';

            busIds = mpcRes.bus(:, 1);
            busVm  = mpcRes.bus(:, 8);
            busVa  = mpcRes.bus(:, 9);
            for k = 1:nIC
                toB = brComplete(idxIC(k), 2);
                bPos = find(busIds == toB, 1);
                intconv_PF.vmod(k)   = busVm(bPos);
                intconv_PF.vangle(k) = busVa(bPos);
                intconv_PF.P(k)      = brComplete(idxIC(k), 16) / Sb;
                intconv_PF.Q(k)      = brComplete(idxIC(k), 17) / Sb;
                intconv_PF.Pdc(k)    = brComplete(idxIC(k), 14) / Sb;
            end

            writeBr = isfield(session.config, 'powerflow') && ...
                      isfield(session.config.powerflow, 'write_branch') && ...
                      session.config.powerflow.write_branch;
            if writeBr || strcmp(solverKind, 'OPF')
                icNames = session.tables.IntConv.Properties.VariableNames(2:end);
                nCopy = min(size(brComplete, 2), numel(icNames));
                for c = 1:nCopy
                    cn = icNames{c};
                    session.tables.IntConv.(cn) = brComplete(idxIC, c);
                    session.rawTables.IntConv = local_update_raw_col(session.rawTables.IntConv, cn, brComplete(idxIC, c));
                end
            end
        end
    end
end

% Build branchFlows table for UI diagram & static report
if ~isempty(brComplete) && size(brComplete, 2) >= 17
    fbus = brComplete(:, 1);
    tbus = brComplete(:, 2);
    Pf   = brComplete(:, 14);
    Qf   = brComplete(:, 15);
    Pt   = brComplete(:, 16);
    Qt   = brComplete(:, 17);
    Ploss = Pf + Pt;
    rateA = brComplete(:, 6);
    Smax = max(hypot(Pf, Qf), hypot(Pt, Qt));
    loading = zeros(size(fbus));
    hasRate = rateA > 0;
    loading(hasRate) = 100.0 * Smax(hasRate) ./ rateA(hasRate);
    if size(brComplete, 2) >= 26
        isIC = brComplete(:, 26) ~= 0;
    else
        isIC = false(size(fbus));
    end
    isDC = (brComplete(:, 4) == 0) & ~isIC;
    branchFlows = table(fbus, tbus, Pf, Qf, Pt, Qt, Ploss, loading, isDC, isIC, ...
        'VariableNames', {'fbus','tbus','Pf_MW','Qf_MVAr','Pt_MW','Qt_MVAr','Ploss_MW','loading_pct','isDC','isIntConv'});
else
    branchFlows = table(zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), false(0,1), ...
        'VariableNames', {'fbus','tbus','Pf_MW','Qf_MVAr','Pt_MW','Qt_MVAr','Ploss_MW','loading_pct','isDC','isIntConv'});
end

summary = struct( ...
    'solver', solverKind, ...
    'Pgen_MW', sum(mpcRes.gen(:, 2)), ...
    'Qgen_MVAr', sum(mpcRes.gen(:, 3)), ...
    'Pload_MW', sum(mpcRes.bus(:, 3)), ...
    'Qload_MVAr', sum(mpcRes.bus(:, 4)), ...
    'Ploss_MW', sum(branchFlows.Ploss_MW), ...
    'Vm_min', min(mpcRes.bus(:, 8)), ...
    'Vm_max', max(mpcRes.bus(:, 8)));
if isfield(mpcRes, 'f')
    summary.cost = mpcRes.f;
end

pfRes = struct( ...
    'converged', converged, ...
    'kind', solverKind, ...
    'mpc', mpcRes, ...
    'PF_branch_complete', brComplete, ...
    'intconv_PF', intconv_PF, ...
    'branchFlows', branchFlows, ...
    'summary', summary, ...
    'elapsed', elapsed);
end

function TRaw = local_update_raw_col(TRaw, colName, newVals)
if isempty(TRaw) || ~ismember(colName, TRaw.Properties.VariableNames) || height(TRaw) ~= numel(newVals)
    return;
end
col = TRaw.(colName);
if isnumeric(col)
    TRaw.(colName) = double(newVals(:));
else
    sCol = strtrim(string(col));
    isSym = startsWith(sCol, "$");
    for k = 1:numel(sCol)
        if ~isSym(k)
            sCol(k) = string(num2str(newVals(k), '%.12g'));
        end
    end
    TRaw.(colName) = sCol;
end
end
