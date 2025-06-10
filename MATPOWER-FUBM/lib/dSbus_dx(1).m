function [dSbus_dVa, dSbus_dVm, dSbus_dPfsh, dSbus_dQfma, dSbus_dBeqz,...
    dSbus_dBeqv, dSbus_dVtma, dSbus_dQtma, dSbus_dPfdp] = dSbus_dx(Ybus, branch, V, vcart)
%%DSBUS_DX  Calls all functions that compute the partial derivatives of power injection w.r.t. x.
%
%   [DSBUS_DVA, DSBUS_DVM, DSBUS_DPFSH, DSBUS_DQFMA, DSBUS_DBEQZ,...
%        DSBUS_DBEQV, DSBUS_DVTMA, DSBUS_DQTMA, DSBUS_DPFDP] = DSBUS_DX(YBUS, BRANCH, V)
%   Returns sparse matrices containing the partial derivatives w.r.t. x.
%
%   where:
%
%      x = [Va, Vm, Theta_Shifter, ma, Beq]
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSbus_dVa, dSbus_dVm, dSbus_dPfsh, dSbus_dQfma, dSbus_dBeqz,...
%        dSbus_dBeqv, dSbus_dVtma, dSbus_dQtma] = dSbus_dx(Ybus, branch,V);
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
[dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V, vcart);
%% Shift Angle Partials
[dSbus_dPfsh] = dSbus_dsh(branch, V, 1, vcart);
%% ma Partials
[dSbus_dQfma] = dSbus_dma(branch, V, 1, vcart);
[dSbus_dQtma] = dSbus_dma(branch, V, 2, vcart);
[dSbus_dVtma] = dSbus_dma(branch, V, 4, vcart);
%% Beq Partials
[dSbus_dBeqz] = dSbus_dBeq(branch, V, 1, vcart);
[dSbus_dBeqv] = dSbus_dBeq(branch, V, 2, vcart);
%% Shift Angle Partials for Voltage Droop Control
[dSbus_dPfdp] = dSbus_dsh(branch, V, 3, vcart);



