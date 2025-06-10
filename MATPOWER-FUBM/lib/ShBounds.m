function [ShAngmin, ShAngmax]  = ShBounds(branch, iPxsh, nPxsh)
%SHBOUNDS  Construct shift angle var limits for elements with active power flow control
%   [SHANGMIN, SHANGMAX]  = SHBOUNDS(BRANCH, IPXSH, NPXSH)
%
%   Constructs the var limits for the shifter angles. iPxsh is the "id"
%   branch row where a shifter is located. nPxsh is the number of Shifters. 
%
%       ShAngmin <= ShAng <= ShAngmax
%
%   The limits are given in the SH_MIN and SH_MAX columns of the branch
%   matrix. Shifter angles are taken to be unbounded below if
%   SH_MIN <= -360 and unbounded above if SH_MAX >= 360. 
%
%   Example:
%       [ShAngmin, ShAngmax]  = ShBounds(branch, iPfsh, nPfsh);

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
ShAngmax =  Inf(nPxsh, 1); %Shifter angle upper limit  INF
ShAngmin = -ShAngmax;      %Shifter angle lower limit -INF

%% Find which Shifters have angle limit
iShl = find((branch(iPxsh, SH_MIN) > -360) & ...
            (branch(iPxsh, SH_MIN) <  360) ); %find which Shifters have lower angle boundary, index Shifter Angle lower
iShu = find((branch(iPxsh, SH_MAX) > -360) & ...
            (branch(iPxsh, SH_MAX) <  360) ); %find which Shifters have upper angle boundary, index Shifter Angle upper
if iShl %if we have Shifters with lower angle control
    ShAngmin(iShl,1) = branch(iPxsh(iShl), SH_MIN)* (pi/180); %Sets lower boundary in radians
end
if iShu %if we have Shifters with upper angle control
    ShAngmax(iShu,1) = branch(iPxsh(iShu), SH_MAX)* (pi/180); %Sets upper boundary in radians
end
