function opfRes = vflexp_run_opf(session)
%VFLEXP_RUN_OPF Execute Optimal Power Flow (AC or hybrid AC/DC FUBM) via MATPOWER.
%
%   opfRes = vflexp_run_opf(session)

tStart = tic;
session.evalTables();

% Build MATPOWER struct
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
onGenBuses = mpc.gen(mpc.gen(:, 8) > 0, 1);
if ~any(~isDCBus(:) & mpc.bus(:, 2) == 3 & ismember(mpc.bus(:, 1), onGenBuses))
    candRef = find(~isDCBus(:) & ismember(mpc.bus(:, 1), onGenBuses), 1);
    if ~isempty(candRef), mpc.bus(candRef, 2) = 3; end
end

% Ensure non-degenerate Pmax/Pmin/Qmax/Qmin bounds if user left 0s or pu values in Excel
for g = 1:size(mpc.gen, 1)
    pg = mpc.gen(g, 2);
    qg = mpc.gen(g, 3);
    if mpc.gen(g, 9) <= mpc.gen(g, 10) || mpc.gen(g, 9) < pg
        mpc.gen(g, 9)  = max([mpc.gen(g, 9), abs(pg) * 1.5, mpc.baseMVA]);
        mpc.gen(g, 10) = min([mpc.gen(g, 10), pg, 0]);
    end
    if mpc.gen(g, 4) <= mpc.gen(g, 5) || mpc.gen(g, 4) < abs(qg)
        qLim = max([abs(mpc.gen(g, 4)), abs(qg) * 1.5, 0.5 * mpc.baseMVA]);
        mpc.gen(g, 4) =  qLim;
        mpc.gen(g, 5) = -qLim;
    end
end
for b = 1:size(mpc.bus, 1)
    if mpc.bus(b, 12) <= mpc.bus(b, 13) || mpc.bus(b, 12) == 0
        mpc.bus(b, 12) = 1.10;
        mpc.bus(b, 13) = 0.90;
    end
end

mpc.gencost = repmat([2, 0, 0, 3, 0.01, 20, 0], height(TGen), 1);

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
        % Avoid narrow 0-degree angle bounds blocking OPF
        badAng = (brLines(:, 12) == 0 & brLines(:, 13) == 0);
        brLines(badAng, 12) = -360;
        brLines(badAng, 13) =  360;
    else
        brLines = zeros(0, 13);
    end

    if hasIC
        brIC = table2array(TIC(:, 2:width(TIC)));
        nFubmCols = size(brIC, 2);
        if hasLines
            padLines = zeros(size(brLines, 1), nFubmCols);
            padLines(:, 1:13) = brLines;
            icNames = TIC.Properties.VariableNames(2:end);
            for pair = {{'TAP_MAX',1},{'TAP_MIN',1},{'K2',1},{'SH_MIN',-360},{'SH_MAX',360}}
                idxC = find(strcmp(icNames, pair{1}{1}), 1);
                if ~isempty(idxC), padLines(:, idxC) = pair{1}{2}; end
            end
            mpc.branch = [padLines; brIC];
        else
            mpc.branch = brIC;
        end
    else
        mpc.branch = brLines;
    end
end

mpopt = mpoption('verbose', 0, 'out.all', 0);
mpcRes = runopf(mpc, mpopt);
elapsed = toc(tStart);

Sb = mpc.baseMVA;
converged = isfield(mpcRes, 'success') && mpcRes.success == 1;
brComplete = mpcRes.branch;
intconv_PF = struct('include', false, 'vmod', [], 'vangle', [], 'P', [], 'Q', [], 'Pdc', [], 'from', [], 'to', []);

if converged
    session.tables.Bus.Vm = mpcRes.bus(:, 8);
    session.tables.Bus.Va = mpcRes.bus(:, 9);
    session.rawTables.Bus = local_update_raw(session.rawTables.Bus, 'Vm', mpcRes.bus(:, 8));
    session.rawTables.Bus = local_update_raw(session.rawTables.Bus, 'Va', mpcRes.bus(:, 9));

    session.tables.Gen.Pg = mpcRes.gen(:, 2);
    session.tables.Gen.Qg = mpcRes.gen(:, 3);
    session.tables.Gen.Vg = mpcRes.gen(:, 6);
    session.rawTables.Gen = local_update_raw(session.rawTables.Gen, 'Pg', mpcRes.gen(:, 2));
    session.rawTables.Gen = local_update_raw(session.rawTables.Gen, 'Qg', mpcRes.gen(:, 3));
    session.rawTables.Gen = local_update_raw(session.rawTables.Gen, 'Vg', mpcRes.gen(:, 6));

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
            for k = 1:nIC
                bPos = find(busIds == brComplete(idxIC(k), 2), 1);
                intconv_PF.vmod(k)   = mpcRes.bus(bPos, 8);
                intconv_PF.vangle(k) = mpcRes.bus(bPos, 9);
                intconv_PF.P(k)      = brComplete(idxIC(k), 16) / Sb;
                intconv_PF.Q(k)      = brComplete(idxIC(k), 17) / Sb;
                intconv_PF.Pdc(k)    = brComplete(idxIC(k), 14) / Sb;
            end
        end
    end
end

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

costVal = NaN;
if isfield(mpcRes, 'f'), costVal = mpcRes.f; end
summary = struct( ...
    'solver', 'OPF', ...
    'Pgen_MW', sum(mpcRes.gen(:, 2)), ...
    'Qgen_MVAr', sum(mpcRes.gen(:, 3)), ...
    'Pload_MW', sum(mpcRes.bus(:, 3)), ...
    'Qload_MVAr', sum(mpcRes.bus(:, 4)), ...
    'Ploss_MW', sum(branchFlows.Ploss_MW), ...
    'Vm_min', min(mpcRes.bus(:, 8)), ...
    'Vm_max', max(mpcRes.bus(:, 8)), ...
    'cost', costVal);

opfRes = struct( ...
    'converged', converged, ...
    'kind', 'OPF', ...
    'mpc', mpcRes, ...
    'PF_branch_complete', brComplete, ...
    'intconv_PF', intconv_PF, ...
    'branchFlows', branchFlows, ...
    'summary', summary, ...
    'elapsed', elapsed);

session.pfResult = opfRes;
session.invalidateDynamicResults();

if converged
    session.log('OPF converged in %.3f s (Cost = %.2f $/h, Pgen = %.2f MW, Ploss = %.3f MW).', ...
        elapsed, costVal, summary.Pgen_MW, summary.Ploss_MW);
else
    session.log('WARNING: MATPOWER OPF did NOT converge.');
end
end

function TRaw = local_update_raw(TRaw, colName, newVals)
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
