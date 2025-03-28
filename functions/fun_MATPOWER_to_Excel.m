function fun_MATPOWER_to_Excel(mpc, filename)
% This function converts a MATPOWER case to an Excel file ready to be read
% by VFlexP.
%
% Input arguments:
% mpc: MATPOWER case, as a struct or string indicating .m file name
% filename (optional): the location of Excel file to write in. If omitted,
% the filename is mpc if mpc is a string and is asked if mpc is a
% struct.
%
% Authors: Ýlvaro Moussa-García and Andrés Tomás-Martín
% Instituto de Investigación Tecnológica, ETSI ICAI
% Universidad Pontificia Comillas

if ischar(mpc)
    mpc = string(mpc);
end

if isstring(mpc) % .m file with MATPOWER case
    if nargin == 1 % only one input argument
        filename = mpc;
    end
    mpc = feval(mpc);
else
    if nargin == 1 % only one input argument
        filename = uiputfile(".xlsx", "Excel file to write");
        if ~isstring(filename)
            return
        end
    end
end

filename = char(filename); % convert string to char

% names of variables
var_names_buses = ["bus_i" "type"	"Pd"	"Qd"	"Gs"	"Bs" ...
    "area"	"Vm"	"Va"	"baseKV"	"zone"	"Vmax"	"Vmin"];
var_names_lines = ["fbus"	"tbus"	"r"	"x"	"b"	"rateA"	"rateB"	...
    "rateC"	"ratio"	"angle"	"status"	"angmin"	"angmax"];
var_names_gen = ["bus"	"Pg"	"Qg"	"Qmax"	"Qmin"	"Vg" ...
    "mBase"	"status"	"Pmax"	"Pmin"	"Pc1"	"Pc2"	"Qc1min" ...
    "Qc1max"	"Qc2min"	"Qc2max"	"ramp_agc"	"ramp_10" ...
    "ramp_30"	"ramp_q"	"apf"];
var_names_intconv = ["fbus"	"tbus"	"r"	"x"	"b"	"rateA"	"rateB"	...
    "rateC"	"ratio"	"angle"	"status" "angmin" "angmax" "PF" "QF" "PT" ...
    "QT" "MU_SF" "MU_ST" "MU_ANGMIN" "MU_ANGMAX" "VF_SET" "VT_SET" ...
    "TAP_MAX" "TAP_MIN" "CONV" "BEQ" "K2" "BEQ_MIN" "BEQ_MAX" "SH_MIN" ...
    "SH_MAX" "GSW" "ALPHA1" "ALPHA2" "ALPHA3" "Column37"];

% Line table
Tline=array2table(mpc.branch(:,1:13),"VariableNames",var_names_lines);
Tline.Complexity=strings(height(Tline),1);
Tline.Complexity(:)= 'static omega';
Tline.Complexity_LF=Tline.Complexity;

% Bus table
Tbus=array2table(mpc.bus(:,1:13),"VariableNames",var_names_buses);
Tbus.Complexity_bus=strings(height(Tbus),1);
Tbus.Complexity_bus(:)= 'static omega';
Tbus.Complexity_load=Tbus.Complexity_bus;
Tbus.Complexity_bus_LF=Tbus.Complexity_bus;
Tbus.Complexity_load_LF=Tbus.Complexity_bus;

% Check hybrid AC/DC (FUBM) formulation
if size(mpc.branch, 2) > 13 % FUBM formulation
    idx_lines = mpc.branch(:,26) == 0;
    idx_lines_AC = mpc.branch(idx_lines,4) ~= 0;
    idx_lines_DC = mpc.branch(idx_lines,4) == 0;
    idx_branch_DC =  mpc.branch(:,4) == 0 & idx_lines;
    DC_buses = unique([mpc.branch(idx_branch_DC,1); ...
        mpc.branch(idx_branch_DC,2)]);
    idx_buses_DC = ismember(Tbus.bus_i, DC_buses);
    Tbus.Complexity_bus(idx_buses_DC) = 'DC static omega';
    Tbus.Complexity_bus_LF(idx_buses_DC) = 'DC static omega';
    Tline = array2table(mpc.branch(idx_lines,(1:13)), ...
        "VariableNames", var_names_lines);
    Tline.Complexity = strings(height(Tline),1);
    Tline.Complexity_PF = strings(height(Tline),1);
    Tline.Complexity(idx_lines_DC) = 'DC static omega';
    Tline.Complexity(idx_lines_AC) = 'static omega';
    Tline.Complexity_PF(idx_lines_DC) = 'DC static omega';
    Tline.Complexity_PF(idx_lines_AC) = 'static omega';
else
    Tline = array2table(mpc.branch, "VariableNames", var_names_lines);
    Tline.Complexity = strings(height(Tline),1);
    Tline.Complexity(:) = 'static omega';
    Tline.Complexity_LF = Tline.Complexity;
end

% Gen table
Tgen=array2table(mpc.gen(:,1:21),"VariableNames",var_names_gen);
Tgen.type=strings(height(Tgen),1);
Tgen.type(:) = 'GFr';
Tgen = movevars(Tgen, 'type', 'After', 'bus');

% Interface converters table (only for AC/DC cases)
if size(mpc.branch, 2) > 13 % FUBM formulation
    idx_intconv = mpc.branch(:,26) ~= 0;
    Nintconv = sum(idx_intconv);
    if Nintconv > 0
        Tintconv = [array2table(mpc.branch(idx_intconv,:), ...
            "VariableNames", var_names_intconv)];
        Tintconv.type = strings(Nintconv, 1);
        Tintconv.type(:) = "DC_GFl";
        Tintconv = movevars(Tintconv, "type", "Before", "fbus");
    else
        Tintconv = array2table(ones(0, length(var_names_intconv)), ...
            "VariableNames", var_names_intconv);
    end
else
    Tintconv = array2table(ones(0, length(var_names_intconv)), ...
        "VariableNames", var_names_intconv);
end

% System parameters table
areas = unique(Tbus.area(~startsWith(Tbus.Complexity_bus, 'DC')));
Fb = areas;
Fb(:) = 50;
Tsysparam=array2table([[mpc.baseMVA; zeros(length(areas)-1, 1)], areas, Fb], "VariableNames", [ "Sb" "Area" "Fb"]);
Tsysparam=Tsysparam(:, ["Sb" "Area" "Fb"]);
writetable(Tsysparam, [filename '.xlsx'], "Sheet", "System parameters")

% Parameter variation table
Tparam=array2table([" " , " ", " ", " ", " "], "VariableNames", ...
    [ "Parameter" "Default_value" "From"  "To" "N_steps"]);
Tparam=Tparam(:, ["Parameter" "Default_value" "From"  "To" "N_steps"]);
writetable(Tparam, [filename '.xlsx'], "Sheet", "Parameter variation")

% Secondary control table
Tsec=array2table(["No SC", 0], "VariableNames", ...
    [ "SC_structure" "Define_SC_in_MATLAB"]);
Tsec=Tsec(:, ["SC_structure" "Define_SC_in_MATLAB"]);

% Tertiary control table
Tter=array2table("No TC", "VariableNames",  "TC_structure");

% write tables in Excel file
writetable(Tgen, [filename '.xlsx'],"Sheet","Gen data")
writetable(Tbus, [filename '.xlsx'],"Sheet","Bus data")
writetable(Tline, [filename '.xlsx'],"Sheet","Line data")
writetable(Tintconv, [filename '.xlsx'],"Sheet","IntConv data")
writetable(Tsysparam, [filename '.xlsx'], "Sheet", "System parameters")
writetable(Tparam, [filename '.xlsx'], "Sheet", "Parameter variation")
writetable(Tsec, [filename '.xlsx'], "Sheet", "Secondary control")
writetable(Tter, [filename '.xlsx'], "Sheet", "Tertiary control")
