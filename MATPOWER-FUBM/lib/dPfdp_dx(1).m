function [dPfdp_dVa, dPfdp_dVm, dPfdp_dPfsh, dPfdp_dQfma, dPfdp_dBeqz,...
    dPfdp_dBeqv, dPfdp_dVtma, dPfdp_dQtma, dPfdp_dPfdp] = dPfdp_dx(branch, Yf, Yt, V, flow, vcart)
%%DPFDP_DX  Calls all functions that compute the partial derivatives of Pf-Vf Droop Control w.r.t. x, eq. Pf - Pfset = Kdp*(Vmf - Vmfset).
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

%% default input args
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

%% constants
nb = length(V);             %% number of buses
nl = size(branch, 1);       %% number of lines

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

%% Derivatives of Voltage Magnitude w.r.t. Voltage magnitude
dVmf_dVm = sparse(zeros(nl,nb));                          % Initialize for speed [nl,nb]
fdp = branch(iPfdp, F_BUS);                               % List of "from" buses with Voltage Droop Control [nPfdp, 1]
Cfdp = sparse(1:nPfdp, fdp, ones(nPfdp, 1), nPfdp, nb);   % connection matrix for line & from buses with Voltage Droop Control [nPfdp, nb]
dVmf_dVm(iPfdp,:)=Cfdp;                                   % Fill derivatives [nl, nb]

%% Power Injection Partials
%%Va, Vm Partials
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
%%Shift Angle Partials
[dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1);
%%ma Partials
[dSf_dQfma, dSt_dQfma] = dSbr_dma(branch, V, 1);
[dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2);
[dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4);
%%Beq Partials
[dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, flow);
[dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2);
%%Shift Angle Partials for Voltage Droop Control
[dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch, V, 3);

%% Voltage Droop Control Partials
%%Partials of Pfdp w.r.t. Va 
dPfdp_dVa = -real(dSf_dVa);
%%Partials of Pfdp w.r.t. Vm 
dPfdp_dVm = -real(dSf_dVm) + Kdp.*( dVmf_dVm );
%%Partials of Pfdp w.r.t. ThetaSh for PST, VSCI and VSCII 
dPfdp_dPfsh = -real(dSf_dPfsh);
%%Partials of Pfdp w.r.t. ma 
dPfdp_dQfma = -real(dSf_dQfma);
dPfdp_dQtma = -real(dSf_dQtma);
dPfdp_dVtma = -real(dSf_dVtma);
%%Partials of Pfdp w.r.t. Beq 
dPfdp_dBeqz = -real(dSf_dBeqz);
dPfdp_dBeqv = -real(dSf_dBeqv);
%%Partials of Pfdp w.r.t. ThetaSh for VSCIII 
dPfdp_dPfdp = -real(dSf_dPfdp);


