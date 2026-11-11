function fld = vflexp_state_field(relPath)
%VFLEXP_STATE_FIELD Convert a relative Simulink Integrator path into a deterministic struct field name.
%   fld = VFLEXP_STATE_FIELD(relPath) maps e.g. 'Electromagnetic/order 8/Int1'
%   to 'Electromagnetic_order_8_Int1', ensuring <= 63 chars with zero collisions.

s = char(string(relPath));
s = strrep(s, sprintf('\n'), '_');
s = strrep(s, sprintf('\r'), '_');
s = regexprep(s, '[^A-Za-z0-9_]', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');

if isempty(s)
    s = 'x0';
elseif ~isletter(s(1))
    s = ['x_' s];
end

if numel(s) > 55
    % Deterministic FNV-1a 32-bit hash suffix to avoid collisions on truncation
    bytes = uint32( unicode2native(char(string(relPath)), 'UTF-8') );
    h = uint32(2166136261);
    for k = 1:numel(bytes)
        h = bitxor(h, bytes(k));
        h = uint32(mod(uint64(h) * uint64(16777619), 4294967296));
    end
    fld = sprintf('%s_%08x', s(1:54), h);
else
    fld = s;
end
end
