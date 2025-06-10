function [dAf_dXxma, dAt_dXxma] = ...
                        dAbr_dma(dFf_dXxma, dFt_dXxma, Ff, Ft)
%DABR_DMA  Partial derivatives of squared flow magnitudes w.r.t ma/tap.
%        [DAF_DXXMA, DAT_DXXMA] = ...
%                       DABR_DMA(DFF_DXXMA, DFT_DXXMA, FF, FT)
%   returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   ma/tap
%
%   The following explains the expressions used to form the matrices:
%
%   Let Af refer to the square of the apparent power at the "from" end of
%   each branch,
%
%       Af = abs(Sf).^2           %     At = abs(St).^2
%          = Sf .* conj(Sf)       %        = St .* conj(St)
%          = Pf.^2 + Qf.^2        %        = Pt.^2 + Qt.^2
%
%   then ...
%
%   Partial w.r.t real power,
%       dAf/dPf = 2 * diag(Pf)    %     dAt/dPt = 2 * diag(Pt)
%
%   Partial w.r.t reactive power,
%       dAf/dQf = 2 * diag(Qf)    %     dAt/dQt = 2 * diag(Qt)
%
%   Partial w.r.t ma 
%       dAf/dma = dAf/dPf * dPf/dma + dAf/dQf * dQf/dma
%       dAt/dma = dAt/dPt * dPt/dma + dAt/dQt * dQt/dma
%
%
%   Examples:
%       %% squared current magnitude
%       not coded yet
%
%       %% squared apparent power flow
%       [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = ...
%               dSbr_dV(branch(il,:), Yf, Yt, V);
%       [dFf_dXxma, dFt_dXxma] = ...
%               dSbr_dma(branch(il,:), V, ctrl, vcart);
%       [df_dXxma, dt_dXxma] = ...
%                        dAbr_dma(dFf_dXxma, dFt_dXxma, Ff, Ft)
%
%       %% squared real power flow
%       not coded yet
%
%   See also DIBR_DV, DSBR_DV, DABR_DV.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%                                           TN4-OPF-Derivatives-Cartesian.pdf

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
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% dimensions
nl = length(Ff);

%%----- partials w.r.t. real and imaginary flows -----
dAf_dFfr = sparse(1:nl, 1:nl, 2 * real(Ff), nl, nl); %dAf/dPf
dAf_dFfi = sparse(1:nl, 1:nl, 2 * imag(Ff), nl, nl); %dAf/dQf
dAt_dFtr = sparse(1:nl, 1:nl, 2 * real(Ft), nl, nl); %dAt/dPt
dAt_dFti = sparse(1:nl, 1:nl, 2 * imag(Ft), nl, nl); %dAt/dQt

%% partials w.r.t. ma/tap
dAf_dXxma = dAf_dFfr * real(dFf_dXxma) + dAf_dFfi * imag(dFf_dXxma);
dAt_dXxma = dAt_dFtr * real(dFt_dXxma) + dAt_dFti * imag(dFt_dXxma);

