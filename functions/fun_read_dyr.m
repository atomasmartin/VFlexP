function dyr_data = fun_read_dyr(dyr_path)
% FUN_READ_DYR Read a PSS/E-style .dyr dynamic data file for VFlexP.
%
% Supports:
%   - Line comments starting with // or # (outside quotes)
%   - Multi-line records terminated by /
%   - Variable definitions: $name = <number> /
%   - Variable references: $name (positional) or PARAM=$name (named)
%   - Positional parameters for GENROU and GENSAL
%   - Named parameters NAME=VALUE and COMPLEXITY="A=X;B=Y"
%   - Bus keys (<bus>) and branch keys (<fbus>-<tbus>)

if ~isfile(dyr_path)
    error('VFlexP:fun_read_dyr:FileNotFound', 'Dynamic data file not found: %s', char(dyr_path));
end

raw_text = fileread(char(dyr_path));
% Remove UTF-8 BOM if present
if ~isempty(raw_text) && double(raw_text(1)) == 65279
    raw_text = raw_text(2:end);
end

[rec_lines, rec_texts] = local_split_records(raw_text);

dyr_data = struct();
dyr_data.source = char(dyr_path);
dyr_data.definitions = struct();
dyr_data.definition_names = {};
dyr_data.records = struct('model', {}, 'device_id', {}, 'bus', {}, ...
    'fbus', {}, 'tbus', {}, 'key', {}, 'role', {}, 'params', {}, ...
    'references', {}, 'complexity', {}, 'line_no', {});
dyr_data.warnings = {};

for i = 1:numel(rec_texts)
    rtext = strtrim(rec_texts{i});
    lno = rec_lines(i);
    if isempty(rtext)
        continue;
    end
    [is_def, def_name, def_val] = local_parse_definition(rtext, lno);
    if is_def
        if ~isfield(dyr_data.definitions, def_name)
            dyr_data.definition_names{end+1} = def_name; %#ok<AGROW>
        end
        dyr_data.definitions.(def_name) = def_val;
        assignin('base', def_name, def_val);
        continue;
    end
    rec = local_parse_record(rtext, lno, dyr_data.definitions);
    dyr_data.records(end+1) = rec; %#ok<AGROW>
end
end

function [rec_lines, rec_texts] = local_split_records(text)
lines = regexp(text, '\r\n|\n|\r', 'split');
rec_lines = [];
rec_texts = {};
buf = '';
start_line = 1;
quote_char = '';

for lineno = 1:numel(lines)
    line = local_strip_comment(lines{lineno});
    if isempty(strtrim(line)) && isempty(strtrim(buf))
        continue;
    end
    if isempty(strtrim(buf))
        start_line = lineno;
        trimmed_line = strtrim(line);
        if ~isempty(regexp(trimmed_line, '^\$[A-Za-z_][A-Za-z0-9_]*\s*=\s*[^/]+$', 'once'))
            rec_lines(end+1) = start_line; %#ok<AGROW>
            rec_texts{end+1} = trimmed_line; %#ok<AGROW>
            continue;
        end
    end
    for k = 1:length(line)
        ch = line(k);
        if ~isempty(quote_char)
            buf = [buf, ch]; %#ok<AGROW>
            if ch == quote_char
                quote_char = '';
            end
        elseif ch == '''' || ch == '"'
            quote_char = ch;
            buf = [buf, ch]; %#ok<AGROW>
        elseif ch == '/'
            r = strtrim(buf);
            buf = '';
            if ~isempty(r)
                rec_lines(end+1) = start_line; %#ok<AGROW>
                rec_texts{end+1} = r; %#ok<AGROW>
            end
            start_line = lineno;
        else
            buf = [buf, ch]; %#ok<AGROW>
        end
    end
    buf = [buf, ' ']; %#ok<AGROW>
end

if ~isempty(strtrim(buf))
    error('VFlexP:fun_read_dyr:UnterminatedRecord', ...
        'Unterminated .dyr record starting at line %d (every record must end with "/").', start_line);
end
end

function out = local_strip_comment(line)
out = '';
quote_char = '';
k = 1;
n = length(line);
while k <= n
    ch = line(k);
    if ~isempty(quote_char)
        out = [out, ch]; %#ok<AGROW>
        if ch == quote_char
            quote_char = '';
        end
    elseif ch == '''' || ch == '"'
        quote_char = ch;
        out = [out, ch]; %#ok<AGROW>
    elseif ch == '#'
        break;
    elseif ch == '/' && k < n && line(k+1) == '/'
        break;
    else
        out = [out, ch]; %#ok<AGROW>
    end
    k = k + 1;
end
end

function [is_def, name, val] = local_parse_definition(record_text, line_no)
is_def = false;
name = '';
val = NaN;
cleaned = regexprep(strtrim(record_text), '\s+', ' ');
tok = regexp(cleaned, '^\$([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$', 'tokens', 'once');
if isempty(tok)
    return;
end
is_def = true;
name = tok{1};
val = local_as_number(tok{2});
if isempty(val)
    error('VFlexP:fun_read_dyr:InvalidDefinition', ...
        'Line %d: definition of $%s must be a number, got "%s".', line_no, name, tok{2});
end
end

function tokens = local_tokenise(record_text)
tokens = {};
curr = '';
quote_char = '';
for k = 1:length(record_text)
    ch = record_text(k);
    if ~isempty(quote_char)
        if ch == quote_char
            quote_char = '';
        else
            curr = [curr, ch]; %#ok<AGROW>
        end
    elseif ch == '''' || ch == '"'
        quote_char = ch;
    elseif ch == ' ' || ch == ',' || ch == sprintf('\t')
        if ~isempty(curr)
            tokens{end+1} = curr; %#ok<AGROW>
            curr = '';
        end
    else
        curr = [curr, ch]; %#ok<AGROW>
    end
end
if ~isempty(quote_char)
    error('VFlexP:fun_read_dyr:UnbalancedQuote', 'Unbalanced quote in record: %s', record_text);
end
if ~isempty(curr)
    tokens{end+1} = curr;
end
end

function num = local_as_number(str_val)
t = strtrim(char(str_val));
low = lower(regexprep(t, '^\+', ''));
if strcmp(low, 'inf') || strcmp(low, 'infinity')
    num = Inf;
    return;
elseif strcmp(low, '-inf') || strcmp(low, '-infinity')
    num = -Inf;
    return;
end
val = str2double(t);
if ~isnan(val)
    num = val;
else
    num = [];
end
end

function rec = local_parse_record(text, line_no, definitions)
raw_tokens = local_tokenise(text);
if numel(raw_tokens) < 2
    error('VFlexP:fun_read_dyr:ShortRecord', ...
        'Line %d: record needs at least a key and a model name.', line_no);
end

n_tok = numel(raw_tokens);
tokens = cell(1, n_tok);
came_from = cell(1, n_tok);
for i = 1:n_tok
    [tokens{i}, came_from{i}] = local_resolve_token(raw_tokens{i}, definitions, line_no);
end

key_str = tokens{1};
model_name = upper(strtrim(tokens{2}));
rest = tokens(3:end);
rest_from = came_from(3:end);

rec = struct();
rec.model = model_name;
rec.device_id = '1';
rec.bus = [];
rec.fbus = [];
rec.tbus = [];
rec.key = '';
rec.role = local_model_role(model_name);
rec.params = struct();
rec.references = struct();
rec.complexity = struct();
rec.line_no = line_no;

br_tok = regexp(key_str, '^(-?\d+)\s*-\s*(-?\d+)$', 'tokens', 'once');
if ~isempty(br_tok)
    rec.fbus = str2double(br_tok{1});
    rec.tbus = str2double(br_tok{2});
    rec.key = sprintf('%d-%d', rec.fbus, rec.tbus);
else
    bnum = str2double(key_str);
    if isnan(bnum)
        error('VFlexP:fun_read_dyr:InvalidKey', ...
            'Line %d: "%s" is not a bus number or fbus-tbus pair.', line_no, key_str);
    end
    rec.bus = round(bnum);
    rec.key = sprintf('%d', rec.bus);
end

if ~isempty(rest) && isempty(regexp(rest{1}, '^[A-Za-z_][A-Za-z0-9_]*=', 'once'))
    rec.device_id = rest{1};
    rest = rest(2:end);
    rest_from = rest_from(2:end);
end

pos_order = local_positional_order(model_name);
pos_idx = 1;
seen_named = false;

for k = 1:numel(rest)
    tok = rest{k};
    src = rest_from{k};
    nv = regexp(tok, '^([A-Za-z_][A-Za-z0-9_]*)=(.*)$', 'tokens', 'once');
    if ~isempty(nv)
        seen_named = true;
        pname = nv{1};
        pval_str = nv{2};
        pnum = local_as_number(pval_str);
        if ~isempty(pnum)
            rec.params.(pname) = pnum;
        else
            rec.params.(pname) = pval_str;
        end
        if ~isempty(src)
            rec.references.(pname) = src;
        end
    else
        if seen_named
            error('VFlexP:fun_read_dyr:PositionalAfterNamed', ...
                'Line %d: positional parameter "%s" follows a named parameter.', line_no, tok);
        end
        if pos_idx > numel(pos_order)
            error('VFlexP:fun_read_dyr:TooManyPositional', ...
                'Line %d: model %s takes %d positional parameters.', line_no, model_name, numel(pos_order));
        end
        pname = pos_order{pos_idx};
        pnum = local_as_number(tok);
        if isempty(pnum)
            error('VFlexP:fun_read_dyr:NonNumericPositional', ...
                'Line %d: parameter %s = "%s" is not a number.', line_no, pname, tok);
        end
        rec.params.(pname) = pnum;
        if ~isempty(src)
            rec.references.(pname) = src;
        end
        pos_idx = pos_idx + 1;
    end
end

if isfield(rec.params, 'COMPLEXITY')
    rec.complexity = local_parse_complexity(rec.params.COMPLEXITY);
end
end

function [resolved_tok, var_name] = local_resolve_token(tok, definitions, line_no)
var_name = '';
resolved_tok = tok;
eq_pos = find(tok == '=', 1, 'first');
if isempty(eq_pos)
    param_prefix = '';
    rhs = tok;
else
    param_prefix = tok(1:eq_pos-1);
    rhs = tok(eq_pos+1:end);
    if isempty(regexp(param_prefix, '^[A-Za-z_][A-Za-z0-9_]*$', 'once'))
        return;
    end
end
m = regexp(rhs, '^\$([A-Za-z_][A-Za-z0-9_]*)$', 'tokens', 'once');
if isempty(m)
    return;
end
vname = m{1};
if ~isfield(definitions, vname)
    error('VFlexP:fun_read_dyr:UndefinedVariable', ...
        'Line %d: $%s is used before being defined.', line_no, vname);
end
val = definitions.(vname);
if evalin('base', sprintf('exist(''%s'', ''var'')', vname))
    try
        base_val = evalin('base', vname);
        if isnumeric(base_val) && isscalar(base_val)
            val = base_val;
        end
    catch
    end
end
val_str = num2str(val, '%.15g');
if isempty(param_prefix)
    resolved_tok = val_str;
else
    resolved_tok = sprintf('%s=%s', param_prefix, val_str);
end
var_name = vname;
end

function order = local_positional_order(model_name)
switch model_name
    case 'GENROU'
        order = {'Tpdo', 'Tppdo', 'Tpqo', 'Tppqo', 'H', 'D', ...
                 'Xd', 'Xq', 'Xpd', 'Xpq', 'Xppd', 'Xl', 'S1', 'S12'};
    case 'GENSAL'
        order = {'Tpdo', 'Tppdo', 'Tppqo', 'H', 'D', ...
                 'Xd', 'Xq', 'Xpd', 'Xppd', 'Xl', 'S1', 'S12'};
    otherwise
        order = {};
end
end

function role = local_model_role(model_name)
switch model_name
    case {'DCBUS', 'ACBUS', 'DCLINE'}
        role = 'node';
    case {'GENROU', 'GENSAL'}
        role = 'machine';
    case {'ST4B', 'ST1'}
        role = 'avr';
    case 'IEEEG1'
        role = 'gov';
    case 'PSS2B'
        role = 'pss';
    otherwise
        role = 'device';
end
end

function cstruct = local_parse_complexity(raw_str)
cstruct = struct();
parts = strsplit(char(string(raw_str)), ';');
for i = 1:numel(parts)
    p = strtrim(parts{i});
    if isempty(p)
        continue;
    end
    kv = strsplit(p, '=');
    if numel(kv) == 2
        k = matlab.lang.makeValidName(upper(strtrim(kv{1})));
        v = upper(strtrim(kv{2}));
        cstruct.(k) = v;
    else
        cstruct.MODE = upper(p);
    end
end
end
