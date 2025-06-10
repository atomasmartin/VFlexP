function [Sf, St] = SbrFlows(branch, Yf, Yt, V)
%SBR   Computes branch power flows
%
%   [SF, ST] = SBR(BRANCH, YF, YT, V)
%
%   Returns two matrices containing the complex branch power flows at 
%   "from" and "to" ends of each branch
%
%   The following explains the expressions used to form the matrices:
%
%   If = Yf * V;
%   It = Yt * V;
%   Sf = diag(Vf) * conj(If) = diag(conj(If)) * Vf
%   St = diag(Vt) * conj(It) = diag(conj(It)) * Vt
%
%   Examples:
%   [Sf, St] = Sbr(branch, Yf, Yt, V);
%
%   For more details on the derivations behind the formulation code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010. [Online]. Available:
%          https://matpower.org/docs/TN2-OPF-Derivatives.pdf
%          doi: 10.5281/zenodo.3237866
                                          
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% Define named indices into bus, gen, branch matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3] = idx_brch;%<<AAB-extra fields for FUBM

%% Define
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);

%% Compute intermediate values
Yfc = conj(Yf);
Ytc = conj(Yt);
Vc = conj(V);
Ifc = Yfc * Vc;     %% conjugate of "from" current
Itc = Ytc * Vc;     %% conjugate of "to" current
%% Compute Branch Power Flows
Sf = V(f) .* Ifc;
St = V(t) .* Itc;

