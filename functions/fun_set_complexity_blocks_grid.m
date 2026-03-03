function fun_set_complexity_blocks_grid(gcb, Tline, Tbus, TbusDC, TlineDC, TgenDC)

block = gcb;
mdlWks = get_param(bdroot, 'ModelWorkspace');
bus_complexity_string = Tbus.Complexity_bus; % Complexity of the bus models
load_complexity_string = Tbus.Complexity_load(Tbus.Pd ~= 0 | Tbus.Qd ~= 0); % Complexity of the load models
bvarname = [matlab.lang.makeValidName(get_param(gcb, 'Name')) '_'];
set_param(gcb, "bvarname", bvarname)

if isempty(Tline)
    set_param([gcb '/AC Grid/Dynamic/Lines to buses'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/Lines'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/Selector4'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/Selector6'], 'Commented', 'on')
else
    set_param([gcb '/AC Grid/Dynamic/Lines to buses'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/Lines'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/Selector4'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/Selector6'], 'Commented', 'off')
end

bus_complexity = zeros(1, height(Tbus));
for i=1:length(bus_complexity)
    cmplxty = find(strcmp(bus_complexity_string(i), {...
        'dynamic', 'static', 'static omega', 'only R'}));
    if isempty(cmplxty)
        bus_complexity(i) = 3; % default
    else
        bus_complexity(i) = cmplxty;
    end
end
load_complexity = zeros(1, sum(Tbus.Pd ~= 0 | Tbus.Qd ~= 0));
if isempty(load_complexity)
    set_param([gcb '/AC Grid/Dynamic/Loads'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/Loads to buses'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/Selector2'], 'Commented', 'on')
    set_param([gcb '/AC Grid/Dynamic/iP'], 'Commented', 'on')

else
    set_param([gcb '/AC Grid/Dynamic/Loads'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/Loads to buses'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/Selector2'], 'Commented', 'off')
    set_param([gcb '/AC Grid/Dynamic/iP'], 'Commented', 'off')
end

for i=1:length(load_complexity)
    cmplxty = find(strcmp(load_complexity_string(i), {...
        'series RL dynamic', 'series RL static', 'series RL static omega', ...
        'PQ', 'dynamic PQ', 'parallel RL dynamic', 'parallel RL static', ...
        'parallel RL static omega', 'constant current'}));
    if isempty(cmplxty)
        load_complexity(i) = 3; % default
    else
        load_complexity(i) = cmplxty;
    end
end


% ------ LINES -------
if ~isempty(Tline)
    line_complexity_string = Tline.Complexity; % Complexity of the line models

    line_complexity = zeros(1, height(Tline));
    for i=1:length(line_complexity)
        cmplxty = find(strcmp(line_complexity_string(i), {...
            'dynamic', 'static', 'static omega'}));
        if isempty(cmplxty)
            line_complexity(i) = 3; % default
        else
            line_complexity(i) = cmplxty;
        end
    end
else
    line_complexity = [];
end

if all(bus_complexity == 3) && all(load_complexity == 3) && all(line_complexity == 3)
    set_param([block '/AC Grid'], 'LabelModeActiveChoice', 'Static')
    [Tbus, ~, Tline, TbusDC, TlineDC, ~, ~, ~, ~] = fun_read_mask_tables(block);
    bvarname = get_param(block, "bvarname");
    [y, ynoloads] = fun_admitt(block, Tbus, Tline);
    assignin(mdlWks, [bvarname 'admitt_matrix'], y) % Assign to model workspace
    assignin(mdlWks, [bvarname 'iadmitt_matrix'], inv(y)) % Assign to model workspace
    assignin(mdlWks, [bvarname 'admitt_matrix_noloads'], ynoloads) % Assign to model workspace
    assignin(mdlWks, [bvarname 'iadmitt_matrix_noloads'], inv(ynoloads)) % Assign to model workspace
else
    set_param([block '/AC Grid'], 'LabelModeActiveChoice', 'Dynamic')
end

if sum(bus_complexity == 1) == 0
    set_param([block '/AC Grid/Dynamic/Buses/dynamic'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Buses/dynamic'], 'commented', 'off')
end
if sum(bus_complexity == 2) == 0
    set_param([block '/AC Grid/Dynamic/Buses/static'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Buses/static'], 'commented', 'off')
end
if sum(bus_complexity == 3) == 0
    set_param([block '/AC Grid/Dynamic/Buses/static omega'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Buses/static omega'], 'commented', 'off')
end
if sum(bus_complexity == 4) == 0
    set_param([block '/AC Grid/Dynamic/Buses/only R'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Buses/only R'], 'commented', 'off')
end


if sum(load_complexity == 1) == 0
    set_param([block '/AC Grid/Dynamic/Loads/series RL dynamic'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/series RL dynamic'], 'commented', 'off')
end
if sum(load_complexity == 2) == 0
    set_param([block '/AC Grid/Dynamic/Loads/series RL static'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/series RL static'], 'commented', 'off')
end
if sum(load_complexity == 3) == 0
    set_param([block '/AC Grid/Dynamic/Loads/series RL static omega'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/series RL static omega'], 'commented', 'off')
end
if sum(load_complexity == 4) == 0
    set_param([block '/AC Grid/Dynamic/Loads/PQ'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/PQ'], 'commented', 'off')
end
if sum(load_complexity == 5) == 0
    set_param([block '/AC Grid/Dynamic/Loads/dynamic PQ'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/dynamic PQ'], 'commented', 'off')
end
if sum(load_complexity == 6) == 0
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL dynamic'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL dynamic'], 'commented', 'off')
end
if sum(load_complexity == 7) == 0
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL static'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL static'], 'commented', 'off')
end
if sum(load_complexity == 8) == 0
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL static omega'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Loads/parallel RL static omega'], 'commented', 'off')
end
if sum(load_complexity == 9) == 0
    try
        set_param([block '/AC Grid/Dynamic/Loads/constant current'], 'commented', 'off')
    catch
    end
else
    set_param([block '/AC Grid/Dynamic/Loads/constant current'], 'commented', 'off')
end

assignin(mdlWks, [bvarname 'bus_complexity'], bus_complexity) % Assign to model workspace
assignin(mdlWks, [bvarname 'load_complexity'], load_complexity) % Assign to model workspace


% Comment the unused models for buses, lines and loads
if sum(line_complexity == 1) == 0
    set_param([block '/AC Grid/Dynamic/Lines/dynamic'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Lines/dynamic'], 'commented', 'off')
end
if sum(line_complexity == 2) == 0
    set_param([block '/AC Grid/Dynamic/Lines/static'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Lines/static'], 'commented', 'off')
end
if sum(line_complexity == 3) == 0
    set_param([block '/AC Grid/Dynamic/Lines/static omega'], 'commented', 'on')
else
    set_param([block '/AC Grid/Dynamic/Lines/static omega'], 'commented', 'off')
end

assignin(mdlWks, [bvarname 'line_complexity'], line_complexity) % Assign to model workspace


%% DC System

busDC_complexity_string = TbusDC.Complexity_bus; % Complexity of the bus models
loadDC_complexity_string = TbusDC.Complexity_load; % Complexity of the load models

busDC_complexity = zeros(1, height(TbusDC));
for i=1:length(busDC_complexity)
    cmplxty = find(strcmp(busDC_complexity_string(i), {'DC dynamic', 'DC static'}));
    if isempty(cmplxty)
        busDC_complexity(i) = 2; % default
    else
        busDC_complexity(i) = cmplxty;
    end
end
assignin(mdlWks, [bvarname 'busDC_complexity'], busDC_complexity) % Assign to model workspace

loadDC_complexity = zeros(1, sum(TbusDC.Pd ~= 0));

if isempty(loadDC_complexity)
    set_param([gcb '/DC Grid/Dynamic/Loads'], 'Commented', 'on')
    set_param([gcb '/DC Grid/Dynamic/Loads to buses'], 'Commented', 'on')
    set_param([gcb '/DC Grid/Dynamic/Selector5'], 'Commented', 'on')
    set_param([gcb '/DC Grid/Dynamic/One1'], 'Commented', 'on')
    %set_param([gcb '/DC Grid/Dynamic/iP'], 'Commented', 'on')

else
    set_param([gcb '/DC Grid/Dynamic/Loads'], 'Commented', 'off')
    set_param([gcb '/DC Grid/Dynamic/Loads to buses'], 'Commented', 'off')
    set_param([gcb '/DC Grid/Dynamic/Selector5'], 'Commented', 'off')
    set_param([gcb '/DC Grid/Dynamic/One1'], 'Commented', 'off')
    %set_param([gcb '/DC Grid/Dynamic/iP'], 'Commented', 'off')
end

if isempty(TgenDC)
    set_param([gcb '/DC Generators'], 'Commented', 'on')
else
    set_param([gcb '/DC Generators'], 'Commented', 'off')
end

for i=1:length(loadDC_complexity)
    cmplxty = find(strcmp(loadDC_complexity_string(i), {'DC ConstR', 'DC ConstP','DC ConstI','DC DynamicP'}));
    if isempty(cmplxty)
        loadDC_complexity(i) = 2; % Constant P by default
    else
        loadDC_complexity(i) = cmplxty;
    end
end
assignin(mdlWks, [bvarname 'loadDC_complexity'], loadDC_complexity) % Assign to model workspace

if ~isempty(TlineDC)
    lineDC_complexity_string = TlineDC.Complexity; % Complexity of the line models
    lineDC_complexity = zeros(1, height(TlineDC));
    for i=1:height(TlineDC)
        cmplxty = find(strcmp(lineDC_complexity_string(i), {'DC dynamic', 'DC static'}));
        if isempty(cmplxty)
            lineDC_complexity(i) = 2; % static by default
        else
            lineDC_complexity(i) = cmplxty;
        end
    end
else
    lineDC_complexity = [];
end
assignin(mdlWks, [bvarname 'lineDC_complexity'], lineDC_complexity) % Assign to model workspace

if all(busDC_complexity == 2) && all(loadDC_complexity == 2) && all(lineDC_complexity == 2)
    set_param([block '/DC Grid'], 'LabelModeActiveChoice', 'Static')
else
    set_param([block '/DC Grid'], 'LabelModeActiveChoice', 'Dynamic')
end

if sum(busDC_complexity == 1) == 0
    set_param([block '/DC Grid/Dynamic/Buses/dynamic'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Buses/dynamic'], 'commented', 'off')
end
if sum(busDC_complexity == 2) == 0
    set_param([block '/DC Grid/Dynamic/Buses/static'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Buses/static'], 'commented', 'off')
end
if sum(lineDC_complexity == 1) == 0
    set_param([block '/DC Grid/Dynamic/Lines/Dynamic'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Lines/Dynamic'], 'commented', 'off')
end
if sum(lineDC_complexity == 2) == 0
    set_param([block '/DC Grid/Dynamic/Lines/Static'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Lines/Static'], 'commented', 'off')
end
if sum(loadDC_complexity == 1) == 0
    set_param([block '/DC Grid/Dynamic/Loads/ConstR'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Loads/ConstR'], 'commented', 'off')
end
if sum(loadDC_complexity == 2) == 0
    set_param([block '/DC Grid/Dynamic/Loads/ConstP'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Loads/ConstP'], 'commented', 'off')
end
if sum(loadDC_complexity == 3) == 0
    set_param([block '/DC Grid/Dynamic/Loads/ConstI'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Loads/ConstI'], 'commented', 'off')
end
if sum(loadDC_complexity == 4) == 0
    set_param([block '/DC Grid/Dynamic/Loads/DynamicP'], 'commented', 'on')
else
    set_param([block '/DC Grid/Dynamic/Loads/DynamicP'], 'commented', 'off')
end
end