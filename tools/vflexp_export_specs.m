function outFiles = vflexp_export_specs(projRoot)
%VFLEXP_EXPORT_SPECS Export deterministic text netlists for VFlexP.slx and Simulink/lib/*.slx.
%   outFiles = VFLEXP_EXPORT_SPECS() writes sorted, geometry-free .spec.txt
%   mirrors into Simulink/spec/ for clean line-by-line git diff traceability.

if nargin < 1 || isempty(projRoot)
    projRoot = fileparts(fileparts(mfilename('fullpath')));
end

specDir = fullfile(projRoot, 'Simulink', 'spec');
if ~isfolder(specDir)
    mkdir(specDir);
end

libRoot = fullfile(projRoot, 'Simulink', 'lib');
addpath(genpath(libRoot));

slxList = {fullfile(projRoot, 'Simulink', 'VFlexP.slx')};
d = dir(fullfile(libRoot, '**', '*.slx'));
[~, ord] = sort({d.name});
d = d(ord);
for i = 1:numel(d)
    slxList{end+1} = fullfile(d(i).folder, d(i).name); %#ok<AGROW>
end

outFiles = cell(numel(slxList), 1);
for i = 1:numel(slxList)
    slxFile = slxList{i};
    [~, mdlName] = fileparts(slxFile);
    wasLoaded = bdIsLoaded(mdlName);
    if ~wasLoaded
        load_system(slxFile);
    end

    specPath = fullfile(specDir, [mdlName '.spec.txt']);
    export_single_bd(mdlName, slxFile, projRoot, specPath);
    outFiles{i} = specPath;

    if ~wasLoaded
        close_system(mdlName, 0);
    end
end
end

function export_single_bd(mdl, slxFile, projRoot, specPath)
fid = fopen(specPath, 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));

relSlx = strrep(erase(slxFile, [projRoot filesep]), '\', '/');
bdType = get_param(mdl, 'BlockDiagramType');
fprintf(fid, '# VFlexP Deterministic Simulink Netlist Specification\n');
fprintf(fid, 'model: %s\n', mdl);
fprintf(fid, 'file:  %s\n', relSlx);
fprintf(fid, 'type:  %s\n\n', bdType);

% Do NOT follow external library links: each library has its own .spec.txt!
blks = find_system(mdl, 'MatchFilter', @Simulink.match.allVariants, ...
    'LookUnderMasks', 'all', 'FollowLinks', 'off', 'IncludeCommented', 'on');
blks(1) = []; % exclude root diagram
blks = sort(blks);

fprintf(fid, '## Blocks (%d)\n', numel(blks));
keyParams = {'Value', 'Gain', 'Multiplication', 'InitialCondition', 'Inputs', 'Operator', ...
    'GotoTag', 'TagVisibility', 'VariableName', 'Indices', 'IndexMode', 'NumberOfDimensions', ...
    'Time', 'Before', 'After', 'UpperLimit', 'LowerLimit', 'Threshold', ...
    'VariantControlMode', 'VariantControl', 'OpenFcn'};

for i = 1:numel(blks)
    b = blks{i};
    cleanPath = strrep(strrep(b, sprintf('\n'), ' '), sprintf('\r'), '');
    bt  = get_param(b, 'BlockType');
    cmt = get_param(b, 'Commented');
    ls  = get_param(b, 'LinkStatus');
    hdr = sprintf('[%s] (%s)', cleanPath, bt);
    if ~strcmp(cmt, 'off')
        hdr = sprintf('%s commented=%s', hdr, cmt);
    end
    if ismember(ls, {'resolved', 'inactive', 'unresolved'})
        ref = strrep(get_param(b, 'ReferenceBlock'), sprintf('\n'), ' ');
        hdr = sprintf('%s link=%s->%s', hdr, ls, ref);
    end
    fprintf(fid, '%s\n', hdr);

    % Export mathematical parameters
    dp = get_param(b, 'ObjectParameters');
    for p = 1:numel(keyParams)
        pn = keyParams{p};
        if isfield(dp, pn)
            try
                val = get_param(b, pn);
                if ischar(val) || isstring(val)
                    vs = strtrim(strrep(char(val), sprintf('\n'), ' '));
                    if ~isempty(vs)
                        fprintf(fid, '  .%s = %s\n', pn, vs);
                    end
                end
            catch
            end
        end
    end

    % Export mask parameters & initialization if masked (and not a resolved external link)
    mo = Simulink.Mask.get(b);
    if ~isempty(mo)
        if ~strcmp(ls, 'resolved')
            initCode = strtrim(strrep(mo.Initialization, sprintf('\n'), '; '));
            if ~isempty(initCode)
                fprintf(fid, '  .MaskInitialization = %s\n', initCode);
            end
        end
        for k = 1:numel(mo.Parameters)
            mp = mo.Parameters(k);
            pv = strtrim(strrep(char(mp.Value), sprintf('\n'), ' '));
            fprintf(fid, '  .mask.%s (%s) = %s\n', mp.Name, mp.Type, pv);
        end
    end
end

% Export connections per parent subsystem (without crossing resolved links)
fprintf(fid, '\n## Connections\n');
parents = [{mdl}; find_system(mdl, 'MatchFilter', @Simulink.match.allVariants, ...
    'LookUnderMasks', 'all', 'FollowLinks', 'off', 'IncludeCommented', 'on', 'BlockType', 'SubSystem')];
parents = sort(unique(parents));

edges = {};
for p = 1:numel(parents)
    par = parents{p};
    if ~strcmp(par, mdl) && strcmp(get_param(par, 'LinkStatus'), 'resolved')
        continue;
    end
    lns = find_system(par, 'FindAll', 'on', 'MatchFilter', @Simulink.match.allVariants, ...
        'LookUnderMasks', 'all', 'FollowLinks', 'off', 'IncludeCommented', 'on', ...
        'SearchDepth', 1, 'Type', 'line');
    for k = 1:numel(lns)
        sp = get_param(lns(k), 'SrcPortHandle');
        dp = get_param(lns(k), 'DstPortHandle');
        if sp < 0, continue; end
        srcBlk = strrep(get_param(sp, 'Parent'), sprintf('\n'), ' ');
        srcPrt = get_param(sp, 'PortNumber');
        for d = 1:numel(dp)
            if dp(d) > 0
                dstBlk = strrep(get_param(dp(d), 'Parent'), sprintf('\n'), ' ');
                dstPrt = get_param(dp(d), 'PortNumber');
                edges{end+1, 1} = sprintf('%s:%d -> %s:%d', srcBlk, srcPrt, dstBlk, dstPrt); %#ok<AGROW>
            end
        end
    end
end
edges = sort(unique(edges));
for e = 1:numel(edges)
    fprintf(fid, '%s\n', edges{e});
end
end
