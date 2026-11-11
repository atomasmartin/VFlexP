function [block_DC_buses, block_buses, block_DC_toM, block_toM, dparams, block_fref, block_idx] = mask_ic_routing(blk, unitName)
%MASK_IC_ROUTING Resolve interlinking DC/AC converter routing from ModelWorkspace DC & fref structs.

block_DC_buses = 1;
block_buses    = 1;
block_DC_toM   = 1;
block_toM      = 1;
dparams        = struct('Rf', 0.001, 'Lf', 0.01);
block_fref     = 1;
block_idx      = 1;

if nargin < 2 || isempty(unitName)
    unitName = get_param(blk, 'Name');
end
unitName = char(unitName);

try
    mdl = bdroot(blk);
    if strcmp(get_param(mdl, 'BlockDiagramType'), 'library')
        return;
    end
    mw = get_param(mdl, 'ModelWorkspace');
    if isempty(mw)
        return;
    end

    DC = [];
    fref = [];
    isLocal = false;
    if mw.hasVariable('DC')
        DC = mw.getVariable('DC');
        if mw.hasVariable('fref')
            fref = mw.getVariable('fref');
        end
        if mw.hasVariable([unitName '_init_local'])
            isLocal = logical(mw.getVariable([unitName '_init_local']));
        end
    elseif mw.hasVariable('Grid_DC')
        DC = mw.getVariable('Grid_DC');
        if mw.hasVariable('Grid_fref')
            fref = mw.getVariable('Grid_fref');
        end
    end

    if isstruct(DC) && isfield(DC, [unitName 's_buses']) && ~isempty(DC.([unitName 's_buses']))
        rawBuses   = DC.([unitName 's_buses']);
        rawDCBuses = DC.([unitName 's_DC_buses']);
        n = numel(rawBuses);
        if isLocal
            block_buses    = 1:n;
            block_DC_buses = 1:n;
            block_toM      = eye(n);
            block_DC_toM   = eye(n);
        else
            block_buses    = rawBuses;
            block_DC_buses = rawDCBuses;
            block_toM      = DC.([unitName 's_toM']);
            block_DC_toM   = DC.([unitName 's_DC_toM']);
        end
        if isfield(DC, [unitName '_dparams'])
            dparams = DC.([unitName '_dparams']);
        end
        if isstruct(fref) && isfield(fref, unitName)
            block_fref = fref.(unitName);
        end
        if isfield(DC, [unitName 's_agent_idx']) && ~isempty(DC.([unitName 's_agent_idx']))
            block_idx = DC.([unitName 's_agent_idx']);
        else
            block_idx = 1:n;
        end
    end
catch
    block_DC_buses = 1;
    block_buses    = 1;
    block_DC_toM   = 1;
    block_toM      = 1;
    dparams        = struct('Rf', 0.001, 'Lf', 0.01);
    block_fref     = 1;
    block_idx      = 1;
end
end
