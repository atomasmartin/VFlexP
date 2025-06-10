function [Beqmin, Beqmax]  = BeqBounds(branch,iBeq,nBeq)
%BEQBOUND  Construct Beq var limits for the selected elements
%   [BEQMIN, BEQMAX]  = BEQBOUNDS(BRANCH, IBEQ, NBEQ)
%
%   Obtains the var limits for the BEQ optimisation variable
%   IBEQ is a vector of indices of branches that have BEQ as an
%   optimisation variable.
%
%       BEQMIN <= BEQ <= BEQMAX
%
%   The limits are given in the BEQ_MIN and BEQ_MAX columns of the branch
%   matrix.
%
%   Example:
%       [Beqmin, Beqmax]  = BeqBounds(branch,iBeq,nBeq);

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

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
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
%% Setting all as unbounded
Beqmax =  Inf(nBeq, 1); %AAB-Beq1 angle upper limit  INF
Beqmin = -Beqmax;       %AAB-Beq1 angle lower limit -INF
%% Find which Beq have limits
iBeql = find((branch(iBeq, BEQ_MIN)) ~= 0); %AAB- find which Beq have lower boundary, index of Beq with lower bound
iBequ = find((branch(iBeq, BEQ_MIN)) ~= 0); %AAB- find which Beq have upper boundary, index of Beq with upper bound
if iBeql %AAB- If we have Beq with lower boundary
    Beqmin(iBeql,1) = branch(iBeq(iBeql), BEQ_MIN); %AAB- Sets lower boundary in [pu]
end
if iBequ %AAB- If we have Beq with upper boundary
    Beqmax(iBequ,1) = branch(iBeq(iBequ), BEQ_MAX); %AAB- Sets upper boundary in [pu]
end
