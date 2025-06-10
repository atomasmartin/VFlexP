function [dSf_dVa, dSf_dVm, dSf_dPfsh, dSf_dQfma, dSf_dBeqz,...
    dSf_dBeqv, dSf_dVtma, dSf_dQtma, dSf_dPfdp,...
    dSt_dVa, dSt_dVm, dSt_dPfsh, dSt_dQfma, dSt_dBeqz,...
    dSt_dBeqv, dSt_dVtma, dSt_dQtma, dSt_dPfdp] = dSbr_dx(branch, Yf, Yt, V, vcart)
%%DSBR_DX  Calls all functions that compute the partial derivatives of Sf and St w.r.t. x.
%
%   [DSF_DVA, DSF_DVM, DSF_DPFSH, DSF_DQFMA, DSF_DBEQZ,...
%        DSF_DBEQV, DSF_DVTMA, DSF_DQTMA, DSF_DPFDP,...
%        DST_DVA, DST_DVM, DST_DPFSH, DST_DQFMA, DST_DBEQZ,...
%        DST_DBEQV, DST_DVTMA, DST_DQTMA, DST_DPFDP] = DSBR_DX(BRANCH, YF, YT, V, VCART)

%   Returns sparse matrices containing the partial derivatives w.r.t. x.
%
%   where:
%
%      x = [Va, Vm, Theta_Shifter, ma, Beq]
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVa, dSf_dVm, dSf_dPfsh, dSf_dQfma, dSf_dBeqz,...
%    dSf_dBeqv, dSf_dVtma, dSf_dQtma, dSf_dPfdp,...
%    dSt_dVa, dSt_dVm, dSt_dPfsh, dSt_dQfma, dSt_dBeqz,...
%    dSt_dBeqv, dSt_dVtma, dSt_dQtma, dSt_dPfdp] = dSbr_dx(branch, Yf, Yt, V, 0)
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
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
%% Shift Angle Partials
[dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1);
%% ma Partials
[dSf_dQfma, dSt_dQfma] = dSbr_dma(branch, V, 1);
[dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2);
[dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4);
%% Beq Partials
[dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1);
[dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2);
%% Shift Angle Partials for Voltage Droop Control
[dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3);



