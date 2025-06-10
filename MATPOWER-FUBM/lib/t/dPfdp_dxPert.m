function [num_dPfdp_dPfsh, num_dPfdp_dQfma, num_dPfdp_dBeqz,...
    num_dPfdp_dBeqv, num_dPfdp_dVtma, num_dPfdp_dQtma, num_dPfdp_dPfdp] = dPfdp_dxPert(baseMVA, bus, branch, Yf, Yt, V, pert, vcart)
%%DPFDP_DX  Calls all functions that compute the Finite differences partial derivatives of Pf-Vf Droop Control w.r.t. x, eq. Pf - Pfset = Kdp*(Vmf - Vmfset).
%
%  [DPFDP_DVA, DPFDP_DVM, DPFDP_DPFSH, DPFDP_DQFMA, DPFDP_DBEQZ,...
%              DPFDP_DBEQV, DPFDP_DVTMA, DPFDP_DQTMA, DPFDP_DPFDP] = DPFDP_DX(BRANCH, YF, YT, V, FLOW, VCART)
%
%   Returns sparse matrices containing the partial derivatives w.r.t. x.
%
%   where:
%
%      x = [Va, Vm, Theta_Shifter, ma, Beq]
%
%   FLOW indicates if the derivatives of Beqz will be considered for Power
%   Flow or for Optimal Power Flow, due to solvability.
%   
%       FLOW = 1; Power Flow
%       FLOW = 3; Optimal Power Flow
%
%   The following explains the expressions used to form the matrices:
%   The Voltage Droop Control equation is given by:
%
%   Pf - Pfset = Kdp.*(Vmf - Vmfset)   or  -Pf + Pfset + Kdp.*(Vmf - Vmfset) = 0
%
%   Polar coordinates:
%   Partials of Pfdp w.r.t. Va 
%   dPfdp_dVa = -real(dSf_dVa)   + 0 + Kdp.*( 0        - 0 );
%   Partials of Pfdp w.r.t. Vm 
%   dPfdp_dVm = -real(dSf_dVm)   + 0 + Kdp.*( dVmf/dVm - 0 );
%   Partials of Pfdp w.r.t. ThetaSh for PST, VSCI and VSCII 
%   dPfdp_dPfsh = -real(dSf_dPfsh)   + 0 + Kdp.*( 0       - 0 );
%   Partials of Pfdp w.r.t. ma 
%   dPfdp_dQfma = -real(dSf_dQfma)   + 0 + Kdp.*( 0       - 0 );
%   dPfdp_dQtma = -real(dSf_dQtma)   + 0 + Kdp.*( 0       - 0 );
%   dPfdp_dVtma = -real(dSf_dVtma)   + 0 + Kdp.*( 0       - 0 );
%   Partials of Pfdp w.r.t. Beq 
%   dPfdp_dBeqz = -real(dSf_dBeqz) + 0 + Kdp.*( 0       - 0 );
%   dPfdp_dBeqv = -real(dSf_dBeqv) + 0 + Kdp.*( 0       - 0 );
%   Partials of Pfdp w.r.t. ThetaSh for VSCIII 
%   dPfdp_dPfdp = -real(dSf_dPfdp)   + 0 + Kdp.*( 0       - 0 );
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVa, dSf_dVm, dSf_dPfsh, dSf_dQfma, dSf_dBeqz,...
%    dSf_dBeqv, dSf_dVtma, dSf_dQtma,...
%    dSt_dVa, dSt_dVm, dSt_dPfsh, dSt_dQfma, dSt_dBeqz,...
%    dSt_dBeqv, dSt_dVtma, dSt_dQtma] = dSbr_dx(branch, Yf, Yt, V, 0)
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
%   This is part of the Flexible Universal Branch Model (FUBM) for Matpower
%   For more info about the model, email: 
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk 

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<FUBM-extra fields for FUBM


%% Va, Vm Partials
%DONE OUTSIDE THIS FUNCTION %[num_dSf_dVa, num_dSf_dVm, num_dSt_dVa, num_dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
%% Shift Angle Partials
[num_dPfdp_dPfsh] = dPfdp_dshPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dPfdp_dPfdp] = dPfdp_dshPert(baseMVA, bus, branch, V, 3, pert, vcart);
%% ma Partials
[num_dPfdp_dQfma] = dPfdp_dmaPert(baseMVA, bus, branch, V, 1, pert, vcart);
[num_dPfdp_dQtma] = dPfdp_dmaPert(baseMVA, bus, branch, V, 2, pert, vcart);
[num_dPfdp_dVtma] = dPfdp_dmaPert(baseMVA, bus, branch, V, 4, pert, vcart);
%% Beq Partials
[num_dPfdp_dBeqz] = dPfdp_dBeqPert(baseMVA, bus, branch, V, 3, pert, vcart);
[num_dPfdp_dBeqv] = dPfdp_dBeqPert(baseMVA, bus, branch, V, 2, pert, vcart);



