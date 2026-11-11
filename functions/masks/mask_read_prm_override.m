function val = mask_read_prm_override(blk, unitName, fieldName, fallbackVal)
%MASK_READ_PRM_OVERRIDE Read popup/mode override from prm_<unitName> in ModelWorkspace if present.

val = fallbackVal;
try
    mdl = bdroot(blk);
    if strcmp(get_param(mdl, 'BlockDiagramType'), 'library')
        return;
    end
    mw = get_param(mdl, 'ModelWorkspace');
    if isempty(mw)
        return;
    end
    varName = ['prm_' char(unitName)];
    if mw.hasVariable(varName)
        s = mw.getVariable(varName);
        if isstruct(s) && isfield(s, fieldName)
            val = s.(fieldName);
        end
    end
catch
end
end
