function [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch)
%MAKEYBUS   Builds the bus admittance matrix and branch admittance matrices.
%   [YBUS, YF, YT] = MAKEYBUS(MPC)
%   [YBUS, YF, YT] = MAKEYBUS(BASEMVA, BUS, BRANCH)
%   
%   Returns the full bus admittance matrix (i.e. for all buses) and the
%   matrices YF and YT which, when multiplied by a complex voltage vector,
%   yield the vector currents injected into each line from the "from" and
%   "to" buses respectively of each line. Does appropriate conversions to p.u.
%   Inputs can be a MATPOWER case struct or individual BASEMVA, BUS and
%   BRANCH values. Bus numbers must be consecutive beginning at 1
%   (i.e. internal ordering).
%
%   See also MAKEJAC, MAKESBUS, EXT2INT.

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

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM-extra fields for FUBM

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeYbus: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end
%% Identify FUBM formulation
if (size(branch,2) < KDP) 
    fubm = 0; %Its not a fubm formulation
else
    fubm = 1; %Its a fubm formulation
end 

%% for each branch, compute the elements of the branch admittance matrix where
%%
%%      | If |   | Yff  Yft |   | Vf |
%%      |    | = |          | * |    |
%%      | It |   | Ytf  Ytt |   | Vt |
%%
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
Bc = stat .* branch(:, BR_B);                           %% line charging susceptance
ma = ones(nl, 1);                               %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
ma(i) = branch(i, TAP);                         %% assign non-zero tap ratios
ShAng = branch(:, SHIFT)* pi/180;               %% Shift angle (in degrees from branch) and changed to radians
tap = ma .* exp(1j* ShAng);                     %% add phase shifters to the tap

%% Build Yff, Yft, Ytf and Ytt
if fubm %FUBM formulation
    Beq= stat .* branch(:, BEQ);                                %%FUBM- VSC Equivalent Reactor for absorbing or supplying reactive power and zero constraint in DC side   
    Gsw= stat .* branch(:, GSW);                                %%FUBM- VSC Switching losses
    k2 = branch(:, K2);                                         %%FUBM- VSC constant depending of how many levels does the VSC is simulating. Default k2 for branches = 1, Default k2 for VSC = sqrt(3)/2
    
    Ytt = Ys + 1j*Bc/2;
    Yff = Gsw+( (Ytt+1j*Beq) ./ ((k2.^2).*tap .* conj(tap))  ); %%FUBM- FUBM formulation
    Yft = - Ys ./ ( k2.*conj(tap) );                            %%FUBM- FUBM formulation
    Ytf = - Ys ./ ( k2.*tap );                                  %%FUBM- FUBM formulation
else %Matpower Formulation
    Ytt = Ys + 1j*Bc/2;
    Yff =  Ytt ./ (tap .* conj(tap));
    Yft = - Ys ./ conj(tap) ;
    Ytf = - Ys ./ tap ;
end
%% compute shunt admittance
%% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
%% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
%% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
%% i.e. Ysh = Psh + j Qsh, so ...
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA; %% vector of shunt admittances

%% bus indices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses

%% for best performance, choose method based on MATLAB vs Octave and size
if nb < 300 || have_fcn('octave')   %% small case OR running on Octave
    %% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    %% at each branch's "from" bus, and Yt is the same for the "to" bus end
    i = [1:nl 1:nl]';                           %% double set of row indices
    Yf = sparse(i, [f; t], [Yff; Yft], nl, nb);
    Yt = sparse(i, [f; t], [Ytf; Ytt], nl, nb);

    %% build Ybus
    Ybus = sparse([f;f;t;t], [f;t;f;t], [Yff;Yft;Ytf;Ytt], nb, nb) + ... %% branch admittances
            sparse(1:nb, 1:nb, Ysh, nb, nb);        %% shunt admittance
else                                %% large case running on MATLAB
    %% build connection matrices
    Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
    Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses

    %% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    %% at each branch's "from" bus, and Yt is the same for the "to" bus end
    Yf = sparse(1:nl, 1:nl, Yff, nl, nl) * Cf + sparse(1:nl, 1:nl, Yft, nl, nl) * Ct;
    Yt = sparse(1:nl, 1:nl, Ytf, nl, nl) * Cf + sparse(1:nl, 1:nl, Ytt, nl, nl) * Ct;

    %% build Ybus
    Ybus = Cf' * Yf + Ct' * Yt + ...            %% branch admittances
            sparse(1:nb, 1:nb, Ysh, nb, nb);    %% shunt admittance
end
