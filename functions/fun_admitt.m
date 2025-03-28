function [y, ynoloads] = fun_admitt(gcb, Tbus, Tline)

% Name:    
% fun_admitt.m
%
% Description: 
% Admittance matriz
%
% Author:  
% Luis Rouco 
% Instituto de Investigacion Tecnologica (I.I.T.)
% Madrid
%
% Date:
% February 6, 1989
%
% Code adapted to MATPOWER format by Andrés Tomás-Martín
% in Date: November 20, 2023

Sb = evalin("base", get_param(gcb, 'Sb'));
nbus = height(Tbus);
busn = Tbus.bus_i;
aux = zeros(nbus, 1);
if ~isempty(Tline)
    for i=1:nbus
        aux(i) = sum(Tline.b(Tline.fbus == busn(i) | Tline.tbus == busn(i)))/2;
    end
end
Rbus = Sb./Tbus.Gs;
Cbus = Tbus.Bs/Sb + aux;
Ploads = Tbus.Pd(Tbus.Pd ~= 0 | Tbus.Qd ~= 0)/Sb;
Qloads = Tbus.Qd(Tbus.Pd ~= 0 | Tbus.Qd ~= 0)/Sb;
% normalise with operating voltage
R_Loads = Tbus.Vm(Tbus.Pd ~= 0 | Tbus.Qd ~= 0).^2.*(Ploads./(Ploads.^2+Qloads.^2));
L_Loads = Tbus.Vm(Tbus.Pd ~= 0 | Tbus.Qd ~= 0).^2.*(Qloads./(Ploads.^2+Qloads.^2));
Yloads = zeros(nbus, 1);
Yloads(Tbus.Pd ~= 0 | Tbus.Qd ~= 0) = 1./(R_Loads + 1i*L_Loads);

y=zeros(nbus,nbus);

if ~isempty(Tline)
    nbranch = height(Tline);
    tbusn = Tline.fbus;
    zbusn = Tline.tbus;
    tap = Tline.ratio;
    tap(tap == 0) = 1;
    phase = deg2rad(Tline.angle);
    r = Tline.r;
    x = Tline.x;
    
    for j=1:nbranch
        i1=find(tbusn(j)==busn);
        i2=find(zbusn(j)==busn);    
        yij=1/(r(j)+1i*x(j));    
        y(i1,i1) = y(i1,i1) + (yij/(tap(j)^2));
        y(i2,i2) = y(i2,i2) + yij;
        y(i1,i2) = y(i1,i2) - yij/(tap(j)*(cos(phase(j))-1i*sin(phase(j))));
        y(i2,i1) = y(i2,i1) - yij/(tap(j)*(cos(phase(j))+1i*sin(phase(j))));
    end
    
end

for j=1:nbus
  i1=j;
  y(i1,i1)=y(i1,i1)+Yloads(i1)+1/(Rbus(i1))+1i*Cbus(i1);
end
ynoloads = y - diag(Yloads);

end
