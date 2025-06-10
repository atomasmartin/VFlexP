function [stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb)
%GETBRANCHDATA   Obtains the requested data from branch.
%   [STAT, CF, CT, K2, TAP, YS, BC, BEQ] = GETBRANCHDATA(BRANCH, NB)
 
%   Returns the requested data from branch. 

%   Example:

%   [stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb)

%   See also MAKEYBUS, MAKEYBUS_AAB.

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM

%% constants
nl = size(branch, 1);       %% number of lines

%% Identify FUBM formulation
if (size(branch,2) < KDP) 
    fubm = 0; %Its not a fubm formulation
    error('getbranchdata: There is missing data in the branch matrix. FUBM formulation')
else
    fubm = 1; %Its a fubm formulation
end 
%--------------------------------------------------------------------------
%% calculations
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
Bc = stat .* branch(:, BR_B);                           %% line charging susceptance
Beq= stat .* branch(:, BEQ);                            %%VSC Equivalent Reactor for absorbing or supplying reactive power and zero constraint in DC side   
ma = ones(nl, 1);                               %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
ma(i) = branch(i, TAP);                         %% assign non-zero tap ratios
ShAngle = branch(:, SHIFT);                     %%AAB- Shift angle in degrees for PST and VSC
tap = ma .* exp(1j*pi/180 * ShAngle);           %% add phase shifters %%<<AAB- Changes Shift angle to radians
k2 = branch(:, K2);                             %% VSC constant depending of how many levels does the VSC is simulating. Default k2 for branches = 1, Default k2 for VSC = sqrt(3)/2

%% build connection matrices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
