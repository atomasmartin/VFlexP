function [mamin, mamax]  = maBounds(branch, iXxma, nXxma)
%MABOUNDS  Construct ma/tap var limits for the selected elements
%   [MAMIN, MAMAX]  = MABOUNDS(BRANCH, IXXMA, NXXMA)
%
%   Constructs the var limits for the ma/tap optimization variable. iXxma is the "id"
%   branch row where a controlled ma/tap is located. nXxma is the number of ma's/tap's. 
%
%       mamin <= ma <= mamax 
%
%   The limits are given in the TAP_MIN and TAP_MAX columns of the branch
%   matrix. 
%
%   Examples:
%       [Qtmamin, Qtmamax]  = maBounds(branch, iQtma, nQtma);
%       Qtmamin <= maQt <= Qtmamax
%       [Vtmamin, Vtmamax]  = maBounds(branch, iVtma, nVtma);

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based on and created for MATPOWER
%   email snoop_and@hotmail.com for more info.

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3] = idx_brch;%<<AAB-extra fields for FUBM
%% Setting all as unbounded
mamax =  Inf(nXxma, 1); %ma/tap upper limit  INF
mamin = -mamax;         %ma/tap lower limit -INF

%% Asign ma's/tap's limits
imal = find((branch(iXxma, TAP_MIN)) ~= (branch(iXxma, TAP_MAX))); %AAB- %index of elements with lower boundary
imau = imal; %index of elements with upper boundary
if imal %if we have ma's/tap's with lower angle control
    mamin(imal,1) = branch(iXxma(imal), TAP_MIN); %Sets lower boundary [pu]
end
if imau %if we have ma's/tap's with upper angle control
    mamax(imau,1) = branch(iXxma(imau), TAP_MAX); %Sets upper boundary [pu]
end
