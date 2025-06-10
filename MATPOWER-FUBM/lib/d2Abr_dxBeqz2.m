function [H13, H23, H31, H32, H33] = d2Abr_dxBeqz2(d2Sbr_dxBeqz2, dSbr_dVa, dSbr_dVm, dSbr_dBeq, Sbr, V, mu)
%D2ABR_DXBEQZ2  Computes 2nd derivatives of |branch flow|^2 w.r.t. BeqzVa, BeqzVm, BeqzBeqz and contrariwise
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%   (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H13, H23, H31, H32, H33] = D2ABR_DXBEQZ2(D2SBR_DXBEQZ2, DSBR_DVA, DSBR_DVM, DSBR_DBEQ, SBR, V, MU)
%
%   Returns 5 matrices containing the 2nd partial derivatives of Abr
%   where:
%   Abr = abs(Sbr).^2
%       = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H13 = HBzVa = d2Abr_dBeqzVa   = d/dBeqz ( (dAbr/dVa  )'*mu )
%   H23 = HBzVm = d2Abr_dBeqzVm   = d/dBeqz ( (dAbr/dVm  )'*mu )
%   H31 = HVaBz = d2Abr_dVaBeqz   = d/dVa   ( (dAbr/dBeqz)'*mu )
%   H32 = HVmBz = d2Abr_dVmBeqz   = d/dVm   ( (dAbr/dBeqz)'*mu )
%   H33 = HBzBz = d2Abr_dBeqz2    = d/dBeqz ( (dAbr/dBeqz)'*mu )
%   
%   H13 = HBzVa = d2Abr_dBeqzVa   =  (  2*real(d2Sbr_dBeqzVa )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dBeqz ))'*mu 
%   H23 = HBzVm = d2Abr_dBeqzVm   =  (  2*real(d2Sbr_dBeqzVm )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dBeqz ))'*mu
%   H31 = HVaBz = d2Abr_dVaBeqz   =  (  2*real(d2Sbr_dVaBeqz )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dVa   ))'*mu
%   H32 = HVmBz = d2Abr_dVmBeqz   =  (  2*real(d2Sbr_dVmBeqz )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dVm   ))'*mu
%   H33 = HBzBz = d2Abr_dBeqz2    =  (  2*real(d2Sbr_dBeqz2  )*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dBeqz ))'*mu
%
%   Example:
%
%      [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch(il,:));
%      d2Sf_dBeqz2 = @(V, mu)d2Sf_dxBeqz2(branch(il,:), V, mu);
%      [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%      [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1);
%      [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Abr_dxBeqz2(d2Sf_dxBeqz2, dSf_dVa, dSf_dVm, dSf_dBeq, Sf, V, mu)
%      [Ht13, Ht23, Ht31, Ht32, Ht33] = d2Abr_dxBeqz2(d2St_dxBeqz2, dSt_dVa, dSt_dVm, dSt_dBeq, St, V, mu)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2.
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

[F13, F23, F31, F32, F33] = d2Sbr_dxBeqz2(V,conj(Sbr).*mu);

d2Abr_dBeqzVa  = 2 * real( F13 + dSbr_dVa.'  * diagmu * conj(dSbr_dBeq) ); 
d2Abr_dBeqzVm  = 2 * real( F23 + dSbr_dVm.'  * diagmu * conj(dSbr_dBeq) );
d2Abr_dVaBeqz  = 2 * real( F31 + dSbr_dBeq.' * diagmu * conj(dSbr_dVa ) );
d2Abr_dVmBeqz  = 2 * real( F32 + dSbr_dBeq.' * diagmu * conj(dSbr_dVm ) );
d2Abr_dBeqz2   = 2 * real( F33 + dSbr_dBeq.' * diagmu * conj(dSbr_dBeq) );
 
H13 = sparse(d2Abr_dBeqzVa);
H23 = sparse(d2Abr_dBeqzVm);
H31 = sparse(d2Abr_dVaBeqz);
H32 = sparse(d2Abr_dVmBeqz);
H33 = sparse(d2Abr_dBeqz2);
