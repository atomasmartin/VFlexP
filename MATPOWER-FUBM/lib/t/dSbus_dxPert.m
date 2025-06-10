function [num_dSbus_dPfsh, num_dSbus_dQfma, num_dSbus_dBeqz,...
    num_dSbus_dBeqv, num_dSbus_dVtma, num_dSbus_dQtma, num_dSbus_dPfdp] = dSbus_dxPert(baseMVA, bus, branch, V, pert, vcart)
%%DSBUS_DXPERT  Calls all functions that compute the Finite differences partial derivatives of power injection w.r.t. xfubm.
%
%   [NUM_DSBUS_DPFSH, NUM_DSBUS_DQFMA, NUM_DSBUS_DBEQZ,...
%        NUM_DSBUS_DBEQV, NUM_DSBUS_DVTMA, NUM_DSBUS_DQTMA, NUM_DSBUS_DPFDP] = DSBUS_DXPERT(BASEMVA, BUS, BRANCH, V, PERT, VCART)

%   Returns sparse matrices containing the partial derivatives w.r.t. xfubm
%   using the Finite differences method.
%
%   where:
%
%      xfubm = [Theta_Shifter, ma, Beq]
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [num_dSbus_dPfsh, num_dSbus_dQfma, num_dSbus_dBeqz,...
%           num_dSbus_dBeqv, num_dSbus_dVtma, num_dSbus_dQtma] = dSbus_dxPert(baseMVA, bus, branch, V, pert, vcart)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TNX]  A. Alvarez-Bustos, "AC/DC FUBM and their Derivatives using FUBM
%          Complex Matrix Notation" MATPOWER Technical Note x, Month 20XX.
%             http://www.pserc.cornell.edu/matpower/
%                                           TNX-OPF-Derivatives-FUBM.pdf   %%AAB- Technical note to be written
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% Va, Vm Partials
%DONE OUTSIDE THIS FUNCTION %[dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V, vcart);
%% Shift Angle Partials
[num_dSbus_dPfsh] = dSbus_dshPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSbus_dPfdp] = dSbus_dshPert(baseMVA, bus, branch, V, 3, pert, vcart);
%% ma Partials
[num_dSbus_dQfma] = dSbus_dmaPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSbus_dQtma] = dSbus_dmaPert(baseMVA, bus, branch, V, 2, pert, vcart);
[num_dSbus_dVtma] = dSbus_dmaPert(baseMVA, bus, branch, V, 4, pert, vcart);
%% Beq Partials
[num_dSbus_dBeqz] = dSbus_dBeqPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSbus_dBeqv] = dSbus_dBeqPert(baseMVA, bus, branch, V, 2, pert, vcart);




