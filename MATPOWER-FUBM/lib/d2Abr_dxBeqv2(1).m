function [H14, H24, H34, H41, H42, H43, H44] = d2Abr_dxBeqv2(d2Sbr_dxBeqv2, dSbr_dVa, dSbr_dVm, dSbr_dBeqz, dSbr_dBeqv, Sbr, V, mu)
%D2ABR_DXBEQV2  Computes 2nd derivatives of |branch flow|^2 w.r.t. BeqvVa, BeqvVm, BeqvBeqz, BeqvBeqv and contrariwise
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%   (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H14, H24, H34, H41, H42, H43, H44] = D2ABR_DXBEQV2(D2SBR_DXBEQV2, DSBR_DVA, DSBR_DVM, DSBR_DBEQZ, DSBR_DBEQV, SBR, V, MU)
%
%   Returns 7 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H14 = HBvVa = d2Abr_dBeqvVa   = d/dBeqv  ( (dAbr/dVa  )'*mu )
%   H24 = HBvVm = d2Abr_dBeqvVm   = d/dBeqv  ( (dAbr/dVm  )'*mu )
%   H34 = HBvBz = d2Abr_dBeqvBeqz = d/dBeqv  ( (dAbr/dBeqz)'*mu )
%   H41 = HVaBv = d2Abr_dVaBeqv   = d/dVa    ( (dAbr/dBeqv)'*mu )
%   H42 = HVmBv = d2Abr_dVmBeqv   = d/dVm    ( (dAbr/dBeqv)'*mu )
%   H43 = HBzBv = d2Abr_dBeqzBeqv = d/dBeqz  ( (dAbr/dBeqv)'*mu )
%   H44 = HBvBv = d2Abr_dBeqvBeqv = d/dBeqv  ( (dAbr/dBeqv)'*mu )
%   
%   H14 = HBvVa = d2Abr_dBeqvVa   =  (  2*real(d2Sbr_dBeqvVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dBeqv ))'*mu 
%   H24 = HBvVm = d2Abr_dBeqvVm   =  (  2*real(d2Sbr_dBeqvVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dBeqv ))'*mu
%   H34 = HBvBz = d2Abr_dBeqvBeqz =  (  2*real(d2Sbr_dBeqvBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dBeqv ))'*mu
%   H41 = HVaBv = d2Abr_dVaBeqv   =  (  2*real(d2Sbr_dVaBeqv  )*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dVa   ))'*mu
%   H42 = HVmBv = d2Abr_dVmBeqv   =  (  2*real(d2Sbr_dVmBeqv  )*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dVm   ))'*mu
%   H43 = HBzBv = d2Abr_dBeqzBeqv =  (  2*real(d2Sbr_dBeqzBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dBeqz ))'*mu
%   H44 = HBvBv = d2Abr_dBeqvBeqv =  (  2*real(d2Sbr_dBeqvBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dBeqv ))'*mu
%
%   Example:
%
%      [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch(il,:));
%      d2Sf_dBeqv2 = @(V, mu)d2Sf_dxBeqv2(branch(il,:), V, mu);
%      [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%      [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1);
%      [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2);
%      [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Abr_dxBeqv2(d2Sf_dxBeqv2, dSf_dVa, dSf_dVm, dSf_dBeqz, dSf_dBeqv, Sf, V, mu)
%      [Ht14, Ht24, Ht34, Ht41, Ht42, Ht43, Ht44] = d2Abr_dxBeqv2(d2St_dxBeqv2, dSt_dVa, dSt_dVm, dSt_dBeqz, dSt_dBeqv, St, V, mu)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXBEQZ2.
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
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
nl = length(mu);

diagmu = sparse(1:nl, 1:nl, mu, nl, nl);

[F14, F24, F34, F41, F42, F43, F44] = d2Sbr_dxBeqv2(V,conj(Sbr).*mu);

d2Abr_dVaBeqv    = 2 * real( F14 + dSbr_dVa.'   * diagmu * conj(dSbr_dBeqv) ); 
d2Abr_dVmBeqv    = 2 * real( F24 + dSbr_dVm.'   * diagmu * conj(dSbr_dBeqv) );
d2Abr_dBeqzBeqv  = 2 * real( F34 + dSbr_dBeqz.' * diagmu * conj(dSbr_dBeqv) );
d2Abr_dBeqvVa    = 2 * real( F41 + dSbr_dBeqv.' * diagmu * conj(dSbr_dVa  ) );
d2Abr_dBeqvVm    = 2 * real( F42 + dSbr_dBeqv.' * diagmu * conj(dSbr_dVm  ) );
d2Abr_dBeqvBeqz  = 2 * real( F43 + dSbr_dBeqv.' * diagmu * conj(dSbr_dBeqz) );
d2Abr_dBeqv2     = 2 * real( F44 + dSbr_dBeqv.' * diagmu * conj(dSbr_dBeqv) );
 
H14 = sparse(d2Abr_dVaBeqv);
H24 = sparse(d2Abr_dVmBeqv);
H34 = sparse(d2Abr_dBeqzBeqv);
H41 = sparse(d2Abr_dBeqvVa);
H42 = sparse(d2Abr_dBeqvVm);
H43 = sparse(d2Abr_dBeqvBeqz);
H44 = sparse(d2Abr_dBeqv2);
