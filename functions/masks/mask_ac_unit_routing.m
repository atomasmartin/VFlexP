function [block_buses, block_idx, block_toM, block_fref] = mask_ac_unit_routing(blk, unitName)
%MASK_AC_UNIT_ROUTING Resolve AC generator bus/agent routing from ModelWorkspace without self-modification.

block_buses = 1;
block_idx   = 1;
block_toM   = 1;
block_fref  = 1;

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

    % Try unprefixed first (standalone / init harness / clean VFlexP), then Grid_ prefix
    prefixes = {'', 'Grid_'};
    found = false;
    for p = 1:numel(prefixes)
        pfx = prefixes{p};
        busVar = [pfx unitName 's_buses'];
        if mw.hasVariable(busVar)
            rawBuses = mw.getVariable(busVar);
            if ~isempty(rawBuses)
                if isempty(pfx) && mw.hasVariable([unitName '_init_local']) && mw.getVariable([unitName '_init_local'])
                    block_buses = 1:numel(rawBuses);
                else
                    block_buses = rawBuses;
                end
                if mw.hasVariable([pfx unitName '_idx'])
                    block_idx = mw.getVariable([pfx unitName '_idx']);
                else
                    block_idx = 1:numel(block_buses);
                end
                if mw.hasVariable([pfx unitName 's_toM'])
                    block_toM = mw.getVariable([pfx unitName 's_toM']);
                else
                    block_toM = eye(numel(block_buses));
                end
                if mw.hasVariable([pfx 'fref'])
                    fref = mw.getVariable([pfx 'fref']);
                    if isstruct(fref) && isfield(fref, unitName)
                        block_fref = fref.(unitName);
                    end
                end
                found = true;
                break;
            end
        end
    end

    if ~found || isempty(block_buses)
        block_buses = 1;
        block_idx   = 1;
        block_toM   = 1;
        block_fref  = 1;
    end
catch
    block_buses = 1;
    block_idx   = 1;
    block_toM   = 1;
    block_fref  = 1;
end
end
