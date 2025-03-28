function [Tbus, Tgen, Tline, TbusDC, TlineDC, Tintconv, ...
    Tbuscomplete, Tlinecomplete, TFb, TgenDC, Tgencomplete] = fun_read_mask_tables(gcb)

% Bus table
block = gcb;
maskObj = Simulink.Mask.get(block);
str_bus = string(eval(get_param(block, 'bus_data')));
if isempty(str_bus)
    return
end
var_names = string({maskObj.getDialogControl('bus_data').Columns(:).Name});
Tbus = array2table(str_bus, "VariableNames", var_names);
aux = zeros(height(Tbus), width(Tbus)-4);
for i=1:width(Tbus)-4
    for j=1:height(Tbus)
        aux(j,i) = evalin("base", str_bus{j,i});
    end
end
Tbus = [array2table(aux, "VariableNames", var_names(1:end-4)) Tbus(:,end-3:end)];
Tbuscomplete = Tbus;
TbusDC = Tbuscomplete(startsWith(Tbuscomplete.Complexity_bus, "DC"), :);
Tbus = Tbuscomplete(~startsWith(Tbuscomplete.Complexity_bus, "DC"), :);
% Line table
str_line = string(eval(get_param(block, 'line_data')));
if isempty(str_line)
    Tline = [];
    TlineDC = [];
    Tlinecomplete = [];
else
    var_names = string({maskObj.getDialogControl('line_data').Columns(:).Name});
    Tline = array2table(str_line, "VariableNames", var_names);
    aux = zeros(height(Tline), width(Tline)-2);
    for i=1:width(Tline)-2
        for j=1:height(Tline)
            aux(j,i) = evalin("base", str_line{j,i});
        end
    end
    Tline = [array2table(aux, "VariableNames", var_names(1:end-2)) Tline(:,end-1:end)];
    Tlinecomplete = Tline;
    TlineDC = Tlinecomplete(startsWith(Tlinecomplete.Complexity, "DC"), :);
    Tline = Tlinecomplete(~startsWith(Tlinecomplete.Complexity, "DC"), :);
end


% Gen table
str_gen = string(eval(get_param(block, 'gen_data')));
if isempty(str_gen)
    return
end
var_names = string({maskObj.getDialogControl('gen_data').Columns(:).Name});
Tgen = array2table(str_gen, "VariableNames", var_names);
aux = zeros(height(Tgen), width(Tgen)-1);
for i=[1 3:width(Tgen)]
    for j=1:height(Tgen)
        if i == 1
            aux(j,i) = evalin("base", str_gen{j,i});
        else
            aux(j,i-1) = evalin("base", str_gen{j,i});
        end
    end
end
Tgencomplete = [array2table(aux(:,1), "VariableNames", var_names(1)) Tgen(:,2) ...
    array2table(aux(:,2:end), "VariableNames", var_names(3:end))];
TgenDC = Tgencomplete(startsWith(Tgencomplete.type, "DC"), :);
Tgen = Tgencomplete(~startsWith(Tgencomplete.type, "DC"), :);

% Interface converters table
str_intconv = string(eval(get_param(block, 'intconv_data')));
if isempty(str_intconv)
    Tintconv = [];
else
    var_names = string({maskObj.getDialogControl('intconv_data').Columns(:).Name});
    Tintconv = array2table(str_intconv, "VariableNames", var_names);
    aux = zeros(height(Tintconv), width(Tintconv)-1);
    for i=2:width(Tintconv)
        for j=1:height(Tintconv)
            aux(j,i-1) = evalin("base", str_intconv{j,i});
        end
    end
    Tintconv = [Tintconv(:,1) array2table(aux, "VariableNames", var_names(2:end))];
end

% Frequency areas table
str_Fb = string(eval(get_param(block, 'Fb_data')));
if isempty(str_Fb)
    TFb = [];
else
    var_names = string({maskObj.getDialogControl('Fb_data').Columns(:).Name});
    TFb = array2table(str_Fb, "VariableNames", var_names);
    aux = zeros(height(TFb), width(TFb));
    for i=1:width(TFb)
        for j=1:height(TFb)
            aux(j,i) = evalin("base", str_Fb{j,i});
        end
    end
    TFb = array2table(aux, "VariableNames", var_names);
end

end