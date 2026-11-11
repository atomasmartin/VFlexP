function code = mask_parse_mode(kind, val)
%MASK_PARSE_MODE Convert string or numeric mode selector to integer variant code.
%   kind = 'stator' : 8 ('order 8'/'EMT8'), 4 ('order 4'/'RMS4'), 2 ('order 2'/'RMS2'), 1 ('order 8 nomatrix')
%   kind = 'onoff'  : 1 ('Dynamic'), 0 ('Static')
%   kind = 'vcc'    : 1 ('Dynamic'), 2 ('Dynamic with limitation'), 3 ('Dynamic VC and static CC'),
%                     4 ('Dynamic with P limitation'), 0 ('Static')
%   kind = 'rlrc'   : 1 ('Dynamic'), 2 ('Static'), 3 ('Static omega')

if isstring(val) || ischar(val)
    s = lower(strtrim(char(val)));
else
    s = '';
end

switch lower(kind)
    case 'stator'
        if isnumeric(val) || islogical(val)
            v = double(val(1));
            if v == 8 || v == 1
                code = 8;
            elseif v == 4 || v == 2
                code = 4;
            elseif v == 3
                code = 2;
            else
                code = 8;
            end
            return;
        end
        switch s
            case {'order 8', 'emt8', 'emt', '8'}
                code = 8;
            case {'order 4', 'rms4', 'rms', '4'}
                code = 4;
            case {'order 2', 'rms2', '2'}
                code = 2;
            case {'order 8 nomatrix', 'order 1', '1'}
                code = 1;
            otherwise
                code = 8;
        end

    case 'onoff'
        if isnumeric(val) || islogical(val)
            v = double(val(1));
            % Note: in 2-option popup {Dynamic|Static}, index 1 = Dynamic, index 2 = Static, 0 = Static
            code = double(v == 1);
            return;
        end
        if any(strcmp(s, {'dynamic', 'emt', 'on', 'true', '1'}))
            code = 1;
        else
            code = 0;
        end

    case 'vcc'
        if isnumeric(val) || islogical(val)
            v = double(val(1));
            if v >= 0 && v <= 4
                code = v;
            else
                code = 0;
            end
            return;
        end
        switch s
            case {'dynamic'}
                code = 1;
            case {'dynamic with limitation', 'dynamic_limiting'}
                code = 2;
            case {'dynamic vc and static cc', 'dynamic_vc'}
                code = 3;
            case {'dynamic with p limitation', 'dynamic_limiting_p'}
                code = 4;
            otherwise
                code = 0; % Static
        end

    case 'rlrc'
        if isnumeric(val) || islogical(val)
            v = double(val(1));
            if v == 1 || v == 2 || v == 3
                code = v;
            else
                code = 1;
            end
            return;
        end
        switch s
            case {'dynamic', 'emt'}
                code = 1;
            case {'static'}
                code = 2;
            case {'static omega', 'static_omega', 'rms'}
                code = 3;
            otherwise
                code = 1;
        end

    otherwise
        code = 0;
end
end
