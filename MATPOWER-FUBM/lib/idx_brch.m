function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch
%IDX_BRCH   Defines constants for named column indices to branch matrix.
%   Example:
%
%   [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%   TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%   ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    branch(4, BR_STATUS) = 0;              % take branch 4 out of service
%    Ploss = branch(:, PF) + branch(:, PT); % compute real power loss vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-11 must be included in input matrix (in case file)
%    1  F_BUS       f, from bus number
%    2  T_BUS       t, to bus number
%    3  BR_R        r, resistance (p.u.)
%    4  BR_X        x, reactance (p.u.)
%    5  BR_B        b, total line charging susceptance (p.u.)
%    6  RATE_A      rateA, MVA rating A (long term rating)
%    7  RATE_B      rateB, MVA rating B (short term rating)
%    8  RATE_C      rateC, MVA rating C (emergency rating)
%    9  TAP         ratio, transformer off nominal turns ratio
%    10 SHIFT       angle, transformer phase shift angle (degrees)
%    11 BR_STATUS   initial branch status, 1 - in service, 0 - out of service
%    12 ANGMIN      minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%    13 ANGMAX      maximum angle difference, angle(Vf) - angle(Vt) (degrees)
%                   (The voltage angle difference is taken to be unbounded below
%                    if ANGMIN < -360 and unbounded above if ANGMAX > 360.
%                    If both parameters are zero, it is unconstrained.)
%
%   columns 14-17 are added to matrix after power flow or OPF solution
%   they are typically not present in the input matrix
%    14 PF          real power injected into "from" end of branch (MW)
%    15 QF          reactive power injected into "from" end of branch (MVAr)
%    16 PT          real power injected into "to" end of branch (MW)
%    17 QT          reactive power injected into "to" end of branch (MVAr)
%
%   columns 18-21 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    18 MU_SF       Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
%    19 MU_ST       Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
%    20 MU_ANGMIN   Kuhn-Tucker multiplier lower angle difference limit (u/degree)
%    21 MU_ANGMAX   Kuhn-Tucker multiplier upper angle difference limit (u/degree)
%
%    FUBM extra columns
%    22 VF_SET      Tap Controlled Voltage (p.u.) (Voltage "from")
%    23 VT_SET      Tap Controlled Voltage (p.u.) (Voltage "to")
%    24 TAP_MAX     Tap Max (p.u.)
%    25 TAP_MIN     Tap Min (p.u.)
%    26 CONV        Converter tipe 1 (zero constraint), 2(VDC control), 3(Droop Control)
%    27 BEQ         Equivalent Suceptance (p.u.)
%    28 K2          Modulation Constant for VSC, otherwise 1
%    29 BEQ_MIN     Beq min value (p.u.)
%    30 BEQ_MAX     Beq max value (p.u.)
%    31 SH_MIN      Shift Angle min value (degrees)
%    32 SH_MAX      Shift Angle max value (degrees)
%    33 GSW         Initial condition for VSC Power Loss (p.u.) (Ploss = Gsw.*Vf.^2 )
%    34 ALPH1       Alpha constant 1 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
%    34 ALPH2       Alpha constant 2 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
%    34 ALPH3       Alpha constant 3 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
%    37 KDP         Constant for VSC Droop Control  (p.u.)   (   Pf - Pfset = Kdp*(Vf - Vfset)  )

%   See also DEFINE_CONSTANTS.

%   ABRAHAM ALVAREZ BUSTOS
%   This code has been modified to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define the indices
F_BUS       = 1;    %% f, from bus number
T_BUS       = 2;    %% t, to bus number
BR_R        = 3;    %% r, resistance (p.u.)
BR_X        = 4;    %% x, reactance (p.u.)
BR_B        = 5;    %% b, total line charging susceptance (p.u.)
RATE_A      = 6;    %% rateA, MVA rating A (long term rating)
RATE_B      = 7;    %% rateB, MVA rating B (short term rating)
RATE_C      = 8;    %% rateC, MVA rating C (emergency rating)
TAP         = 9;    %% ratio, transformer off nominal turns ratio
SHIFT       = 10;   %% angle, transformer phase shift angle (degrees)
BR_STATUS   = 11;   %% initial branch status, 1 - in service, 0 - out of service
ANGMIN      = 12;   %% minimum angle difference, angle(Vf) - angle(Vt) (degrees)
ANGMAX      = 13;   %% maximum angle difference, angle(Vf) - angle(Vt) (degrees)

%% included in power flow solution, not necessarily in input
PF          = 14;   %% real power injected at "from" bus end (MW)       (not in PTI format)
QF          = 15;   %% reactive power injected at "from" bus end (MVAr) (not in PTI format)
PT          = 16;   %% real power injected at "to" bus end (MW)         (not in PTI format)
QT          = 17;   %% reactive power injected at "to" bus end (MVAr)   (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_SF       = 18;   %% Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 19;   %% Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
MU_ANGMIN   = 20;   %% Kuhn-Tucker multiplier lower angle difference limit (u/degree)
MU_ANGMAX   = 21;   %% Kuhn-Tucker multiplier upper angle difference limit (u/degree)
%% fubm extra columns
VF_SET      = 22;   %% Tap Controlled Voltage (p.u.) (Voltage "from")
VT_SET      = 23;   %% Tap Controlled Voltage (p.u.) (Voltage "to")
TAP_MAX     = 24;   %% Tap Max (p.u.)
TAP_MIN     = 25;   %% Tap Min (p.u.)
CONV        = 26;   %% Converter tipe 1 (zero constraint), 2(VDC control), 3(Droop Control)
BEQ         = 27;   %% Equivalent Suceptance (p.u)
K2          = 28;   %% Modulation Constant for VSC, otherwise 1
BEQ_MIN     = 29;   %% Beq min value (p.u.)
BEQ_MAX     = 30;   %% Beq max value (p.u.)
SH_MIN      = 31;   %% Shift Angle min value (degrees)
SH_MAX      = 32;   %% Shift Angle max value (degrees)
GSW         = 33;   %% Initial condition for VSC Power Loss (p.u.) (Ploss = Gsw.*Vf.^2 )
ALPH1       = 34;   %% Alpha constant 1 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
ALPH2       = 35;   %% Alpha constant 2 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
ALPH3       = 36;   %% Alpha constant 3 for VSC Power Loss (Ploss = ALPH3*It^2 + ALPH2*It + ALPH1)
KDP         = 37;   %% Constant for VSC Droop Control (p.u.)     (   Pf - Pfset = Kdp*(Vf - Vfset)  ) 