function [num_dSf_dPfsh, num_dSf_dQfma, num_dSf_dBeqz,...
    num_dSf_dBeqv, num_dSf_dVtma, num_dSf_dQtma, num_dSf_dPfdp,...
    num_dSt_dPfsh, num_dSt_dQfma, num_dSt_dBeqz,...
    num_dSt_dBeqv, num_dSt_dVtma, num_dSt_dQtma, num_dSt_dPfdp] = dSbr_dxPert(baseMVA, bus, branch, Yf, Yt, V, pert, vcart)
%%DSBR_DXPERT  Calls all functions that compute the Finite differences partial derivatives of Sf and St w.r.t. xfubm.
%
%   [NUM_DSF_DPFSH, NUM_DSF_DQFMA, NUM_DSF_DBEQZ,...
%       NUM_DSF_DBEQV, NUM_DSF_DVTMA, NUM_DSF_DQTMA,...
%       NUM_DST_DPFSH, NUM_DST_DQFMA, NUM_DST_DBEQZ,...
%       NUM_DST_DBEQV, NUM_DST_DVTMA, NUM_DST_DQTMA] = DSBR_DXPERT(BASEMVA, BUS, BRANCH, YF, YT, V, PERT, VCART)

%   Returns sparse matrices containing the partial derivatives w.r.t. xfubm.
%
%   where:
%
%      xfubm = [Theta_Shifter, ma, Beq]
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [num_dSf_dPfsh, num_dSf_dQfma, num_dSf_dBeqz,...
%           num_dSf_dBeqv, num_dSf_dVtma, num_dSf_dQtma,...
%           num_dSt_dPfsh, num_dSt_dQfma, num_dSt_dBeqz,...
%           num_dSt_dBeqv, num_dSt_dVtma, num_dSt_dQtma] = dSbr_dxPert(baseMVA, bus, branch, Yf, Yt, V, pert, vcart)
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
%DONE OUTSIDE THIS FUNCTION %[num_dSf_dVa, num_dSf_dVm, num_dSt_dVa, num_dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
%% Shift Angle Partials
[num_dSf_dPfsh, num_dSt_dPfsh] = dSbr_dshPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSf_dPfdp, num_dSt_dPfdp] = dSbr_dshPert(baseMVA, bus, branch, V, 3, pert, vcart);
%% ma Partials
[num_dSf_dQfma, num_dSt_dQfma] = dSbr_dmaPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSf_dQtma, num_dSt_dQtma] = dSbr_dmaPert(baseMVA, bus, branch, V, 2, pert, vcart);
[num_dSf_dVtma, num_dSt_dVtma] = dSbr_dmaPert(baseMVA, bus, branch, V, 4, pert, vcart);
%% Beq Partials
[num_dSf_dBeqz, num_dSt_dBeqz] = dSbr_dBeqPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dSf_dBeqv, num_dSt_dBeqv] = dSbr_dBeqPert(baseMVA, bus, branch, V, 2, pert, vcart);




