function [H18, H28, H38, H48, H58, H68, H78, H81, H82, H83, H84, H85, H86, H87, H88] = d2Abr_dxshdp2(d2Sbr_dxshdp2, dSbr_dVa, dSbr_dVm, dSbr_dBeqz, dSbr_dBeqv, dSbr_dPfsh, dSbr_dqtma, dSbr_dvtma, dSbr_dshdp, Sbr, V, mu)
%D2ABR_DXSHDP2  Computes 2nd derivatives of |branch flow|^2 w.r.t. shdpVa, shdpVm, shdpBeqz, shdpBeqv, shdpSh, shdpqtma, shdpvtma, Vashdp, Vmshdp, Beqzshdp, Beqvshdp, Shshdp, qtmashdp, vtmashdp, shdpshdp.
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H18, H28, H38, H48, H58, H68, H78, H81, H82, H83, H84, H85, H86, H87, H88] = D2ABR_DXSHDP2(D2SBR_DXSHDP2, DSBR_DVA, DSBR_DVM, DSBR_DBEQZ, DSBR_DBEQV, DSBR_DPFSH, DSBR_DQTMA, DSBR_DVTMA, DSBR_DSHDP, SBR, V, MU)
%
%   Returns 15 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H18 = HshdpVa   = d2Abr_dshdpVa     = d/dvtma ( (dAbr/dVa    )'*mu )
%   H28 = HshdpVm   = d2Abr_dshdpVm     = d/dvtma ( (dAbr/dVm    )'*mu )
%   H38 = HshdpBz   = d2Abr_dshdpBeqz   = d/dvtma ( (dAbr/dBeqz  )'*mu )
%   H48 = HshdpBv   = d2Abr_dshdpBeqv   = d/dvtma ( (dAbr/dBeqv  )'*mu )
%   H58 = HshdpSh   = d2Abr_dshdpSh     = d/dvtma ( (dAbr/dSh    )'*mu )
%   H68 = Hshdpqtma = d2Abr_dshdpqtma   = d/dvtma ( (dAbr/dqtma  )'*mu )
%   H78 = Hshdpvtma = d2Abr_dshdpqtma   = d/dvtma ( (dAbr/dvtma  )'*mu )
%   H81 = HVashdp   = d2Abr_dVashdp     = d/dVa   ( (dAbr/dshdp  )'*mu )
%   H82 = HVmshdp   = d2Abr_dVmshdp     = d/dVm   ( (dAbr/dshdp  )'*mu )
%   H83 = HBzshdp   = d2Abr_dBeqzshdp   = d/dBeqz ( (dAbr/dshdp  )'*mu )
%   H84 = HBvshdp   = d2Abr_dBeqvshdp   = d/dBeqv ( (dAbr/dshdp  )'*mu )
%   H85 = HShshdp   = d2Abr_dShshdp     = d/dsh   ( (dAbr/dshdp  )'*mu )
%   H86 = Hqtmashdp = d2Abr_dqtmashdp   = d/dqtma ( (dAbr/dshdp  )'*mu )
%   H87 = Hvtmashdp = d2Abr_dqtmashdp   = d/dvtma ( (dAbr/dshdp  )'*mu )
%   H88 = Hshdpshdp = d2Abr_dshdpshdp   = d/dshdp ( (dAbr/dshdp  )'*mu )
%   
%   H18 = HshdpVa   = d2Abr_dshdpVa     =  (  2*real(d2Sbr_dshdpVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dshdp   ))'*mu 
%   H28 = HshdpVm   = d2Abr_dshdpVm     =  (  2*real(d2Sbr_dshdpVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dshdp   ))'*mu
%   H38 = HshdpBz   = d2Abr_dshdpBeqz   =  (  2*real(d2Sbr_dshdpBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dshdp   ))'*mu
%   H48 = HshdpBv   = d2Abr_dshdpBeqv   =  (  2*real(d2Sbr_dshdpBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dshdp   ))'*mu
%   H58 = HshdpSh   = d2Abr_dshdpSh     =  (  2*real(d2Sbr_dshdpSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dshdp   ))'*mu
%   H68 = Hshdpqtma = d2Abr_dshdpqtma   =  (  2*real(d2Sbr_dshdpqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dshdp   ))'*mu
%   H78 = Hshdpvtma = d2Abr_dshdpvtma   =  (  2*real(d2Sbr_dshdpvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dshdp   ))'*mu
%   H81 = HVashdp   = d2Abr_dVashdp     =  (  2*real(d2Sbr_dVashdp  )*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dVa     ))'*mu
%   H82 = HVmshdp   = d2Abr_dVmshdp     =  (  2*real(d2Sbr_dVmshdp  )*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dVm     ))'*mu
%   H83 = HBzshdp   = d2Abr_dBeqzshdp   =  (  2*real(d2Sbr_dBeqzshdp)*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dBeqz   ))'*mu
%   H84 = HBvshdp   = d2Abr_dBeqvshdp   =  (  2*real(d2Sbr_dBeqvshdp)*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dBeqv   ))'*mu
%   H85 = HShshdp   = d2Abr_dShshdp     =  (  2*real(d2Sbr_dShshdp  )*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dSh     ))'*mu 
%   H86 = Hqtmashdp = d2Abr_dqtmashdp   =  (  2*real(d2Sbr_dqtmashdp)*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dqtma   ))'*mu 
%   H87 = Hvtmashdp = d2Abr_dvtmashdp   =  (  2*real(d2Sbr_dvtmashdp)*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dvtma   ))'*mu 
%   H88 = Hshdpshdp = d2Abr_dshdpshdp   =  (  2*real(d2Sbr_dshdpshdp)*conj(Sbr)  +  dSbr_dshdp *conj(dSbr_dshdp   ))'*mu 
%
%   Example:
%
%      [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch(il,:));
%      d2Sf_dshdp2 = @(V, mu)d2Sf_dxshdp2(branch(il,:), V, mu, mpopt.opf.v_cartesian);
%      [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%      [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1);
%      [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2);
%      [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch(il,:), V, 1);
%      [dSf_dQtma, dSt_dQtma] = dSbr_dqtma(branch(il,:), V, 2);
%      [dSf_dVtma, dSt_dVtma] = dSbr_dvtma(branch(il,:), V, 4);
%      [dSf_dPfdp, dSt_dPfdp] = dSbr_dsh(branch(il,:), V, 3);
%      [Hf18, Hf28, Hf38, Hf48, Hf58, Hf68, Hf78, Hf81, Hf82, Hf83, Hf84, Hf85, Hf86, Hf87, Hf88] = d2Abr_dxshdp2(d2Sf_dxshdp2, dSf_dVa, dSf_dVm, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dqtma, dSf_dvtma, dSf_dshdp, Sf, V, mu)
%      [Ht18, Ht28, Ht38, Ht48, Ht58, Ht68, Ht78, Ht81, Ht82, Ht83, Ht84, Ht85, Ht86, Ht87, Ht88] = d2Abr_dxshdp2(d2St_dxshdp2, dSt_dVa, dSt_dVm, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dqtma, dSt_dvtma, dSf_dshdp, St, V, mu)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXSH2.
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

[F18, F28, F38, F48, F58, F68, F78, F81, F82, F83, F84, F85, F86, F87, F88] = d2Sbr_dxshdp2(V,conj(Sbr).*mu);

d2Abr_dshdpVa    = 2 * real( F18 + dSbr_dVa.'   * diagmu * conj(dSbr_dshdp) ); 
d2Abr_dshdpVm    = 2 * real( F28 + dSbr_dVm.'   * diagmu * conj(dSbr_dshdp) );
d2Abr_dshdpBeqz  = 2 * real( F38 + dSbr_dBeqz.' * diagmu * conj(dSbr_dshdp) );
d2Abr_dshdpBeqv  = 2 * real( F48 + dSbr_dBeqv.' * diagmu * conj(dSbr_dshdp) );
d2Abr_dshdpSh    = 2 * real( F58 + dSbr_dPfsh.' * diagmu * conj(dSbr_dshdp) );
d2Abr_dshdpqtma  = 2 * real( F68 + dSbr_dqtma.' * diagmu * conj(dSbr_dshdp) );
d2Abr_dshdpvtma  = 2 * real( F78 + dSbr_dvtma.' * diagmu * conj(dSbr_dshdp) );
d2Abr_dVashdp    = 2 * real( F81 + dSbr_dshdp.' * diagmu * conj(dSbr_dVa  ) );
d2Abr_dVmshdp    = 2 * real( F82 + dSbr_dshdp.' * diagmu * conj(dSbr_dVm  ) );
d2Abr_dBeqzshdp  = 2 * real( F83 + dSbr_dshdp.' * diagmu * conj(dSbr_dBeqz) );
d2Abr_dBeqvshdp  = 2 * real( F84 + dSbr_dshdp.' * diagmu * conj(dSbr_dBeqv) );
d2Abr_dShshdp    = 2 * real( F85 + dSbr_dshdp.' * diagmu * conj(dSbr_dPfsh) );
d2Abr_dqtmashdp  = 2 * real( F86 + dSbr_dshdp.' * diagmu * conj(dSbr_dqtma) );
d2Abr_dvtmashdp  = 2 * real( F87 + dSbr_dshdp.' * diagmu * conj(dSbr_dvtma) );
d2Abr_dshdp2     = 2 * real( F88 + dSbr_dshdp.' * diagmu * conj(dSbr_dshdp) );
 
H18 = sparse(d2Abr_dshdpVa);
H28 = sparse(d2Abr_dshdpVm);
H38 = sparse(d2Abr_dshdpBeqz);
H48 = sparse(d2Abr_dshdpBeqv);
H58 = sparse(d2Abr_dshdpSh);
H68 = sparse(d2Abr_dshdpqtma);
H78 = sparse(d2Abr_dshdpvtma);
H81 = sparse(d2Abr_dVashdp);
H82 = sparse(d2Abr_dVmshdp);
H83 = sparse(d2Abr_dBeqzshdp);
H84 = sparse(d2Abr_dBeqvshdp);
H85 = sparse(d2Abr_dShshdp);
H86 = sparse(d2Abr_dqtmashdp);
H87 = sparse(d2Abr_dvtmashdp);
H88 = sparse(d2Abr_dshdp2);
