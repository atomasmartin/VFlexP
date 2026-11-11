function [block_buses, block_Pg, block_toM, block_idx] = mask_dc_unit_routing(blk, unitName)
%MASK_DC_UNIT_ROUTING Resolve DC generator routing from ModelWorkspace DC struct.

block_buses = 1;
block_Pg    = 0;
block_toM   = 1;
block_idx   = 1;

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
    isLocal = false;
    if mw.hasVariable('DC')
        DC = mw.getVariable('DC');
        if mw.hasVariable([unitName '_init_local'])
            isLocal = logical(mw.getVariable([unitName '_init_local']));
        end
    elseif mw.hasVariable('Grid_DC')
        DC = mw.getVariable('Grid_DC');
    end

    if isstruct(DC) && isfield(DC, [unitName 's_buses']) && ~isempty(DC.([unitName 's_buses']))
        rawBuses = DC.([unitName 's_buses']);
        if isLocal
            block_buses = 1:numel(rawBuses);
            block_toM   = eye(numel(rawBuses));
        else
            block_buses = rawBuses;
            if isfield(DC, [unitName 's_toM'])
                block_toM = DC.([unitName 's_toM']);
            else
                block_toM = eye(numel(rawBuses));
            end
        end
        if isfield(DC, [unitName 's_Pg'])
            block_Pg = DC.([unitName 's_Pg']);
        else
            block_Pg = zeros(1, numel(rawBuses));
        end
        if isfield(DC, [unitName 's_agent_idx']) && ~isempty(DC.([unitName 's_agent_idx']))
            block_idx = DC.([unitName 's_agent_idx']);
        else
            block_idx = 1:numel(rawBuses);
        end
    end
catch
    block_buses = 1;
    block_Pg    = 0;
    block_toM   = 1;
    block_idx   = 1;
end
end
