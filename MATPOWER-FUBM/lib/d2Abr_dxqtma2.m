function [H16, H26, H36, H46, H56, H61, H62, H63, H64, H65, H66] = d2Abr_dxqtma2(d2Sbr_dxqtma2, dSbr_dVa, dSbr_dVm, dSbr_dBeqz, dSbr_dBeqv, dSbr_dPfsh, dSbr_dQtma, Sbr, V, mu)
%D2ABR_DXQTMA2  Computes 2nd derivatives of |branch flow|^2 w.r.t. qtmaVa, qtmaVm, qtmaBeqz, qtmaBeqv, qtmaSh, Vaqtma, Vmqtma, Beqzqtma, Beqvqtma, Shqtma, qtmaqtma.
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H16, H26, H36, H46, H56, H61, H62, H63, H64, H65, H66] = D2ABR_DXQTMA2(D2SBR_DXMA2, DSBR_DVA, DSBR_DVM, DSBR_DBEQZ, DSBR_DBEQV, DSBR_DPFSH, DSBR_DQTMA, SBR, V, MU)
%
%   Returns 11 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H16 = HqtmaVa   = d2Abr_dqtmaVa     = d/dqtma ( (dAbr/dVa    )'*mu )
%   H26 = HqtmaVm   = d2Abr_dqtmaVm     = d/dqtma ( (dAbr/dVm    )'*mu )
%   H36 = HqtmaBz   = d2Abr_dqtmaBeqz   = d/dqtma ( (dAbr/dBeqz  )'*mu )
%   H46 = HqtmaBv   = d2Abr_dqtmaBeqv   = d/dqtma ( (dAbr/dBeqv  )'*mu )
%   H56 = HqtmaSh   = d2Abr_dqtmaSh     = d/dqtma ( (dAbr/dSh    )'*mu )
%   H61 = HVaqtma   = d2Abr_dVaqtma     = d/dVa   ( (dAbr/dqtma  )'*mu )
%   H62 = HVmqtma   = d2Abr_dVmqtma     = d/dVm   ( (dAbr/dqtma  )'*mu )
%   H63 = HBzqtma   = d2Abr_dBeqzqtma   = d/dBeqz ( (dAbr/dqtma  )'*mu )
%   H64 = HBvqtma   = d2Abr_dBeqvqtma   = d/dBeqv ( (dAbr/dqtma  )'*mu )
%   H65 = HShqtma   = d2Abr_dShqtma     = d/dsh   ( (dAbr/dqtma  )'*mu )
%   H66 = Hqtmaqtma = d2Abr_dqtmaqtma   = d/dqtma ( (dAbr/dqtma  )'*mu )
%   
%   H16 = HqtmaVa   = d2Abr_dqtmaVa     =  (  2*real(d2Sbr_dqtmaVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dqtma   ))'*mu 
%   H26 = HqtmaVm   = d2Abr_dqtmaVm     =  (  2*real(d2Sbr_dqtmaVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dqtma   ))'*mu
%   H36 = HqtmaBz   = d2Abr_dqtmaBeqz   =  (  2*real(d2Sbr_dqtmaBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dqtma   ))'*mu
%   H46 = HqtmaBv   = d2Abr_dqtmaBeqv   =  (  2*real(d2Sbr_dqtmaBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dqtma   ))'*mu
%   H56 = HqtmaSh   = d2Abr_dqtmaSh     =  (  2*real(d2Sbr_dqtmaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dqtma   ))'*mu
%   H61 = HVaqtma   = d2Abr_dVaqtma     =  (  2*real(d2Sbr_dVaqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dVa     ))'*mu
%   H62 = HVmqtma   = d2Abr_dVmqtma     =  (  2*real(d2Sbr_dVmqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dVm     ))'*mu
%   H63 = HBzqtma   = d2Abr_dBeqzqtma   =  (  2*real(d2Sbr_dBeqzqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dBeqz   ))'*mu
%   H64 = HBvqtma   = d2Abr_dBeqvqtma   =  (  2*real(d2Sbr_dBeqvqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dBeqv   ))'*mu
%   H65 = HShqtma   = d2Abr_dShqtma     =  (  2*real(d2Sbr_dShqtma  )*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dSh     ))'*mu 
%   H66 = Hqtmaqtma = d2Abr_dqtmaqtma   =  (  2*real(d2Sbr_dqtmaqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dma     ))'*mu 
%
%   Example:
%
%      [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch(il,:));
%      d2Sf_dqtma2 = @(V, mu)d2Sf_dxqtma2(branch(il,:), V, mu, mpopt.opf.v_cartesian);
%      [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%      [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1);
%      [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2);
%      [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch(il,:), V, 1);
%      [dSf_dQtma, dSt_dQtma] = dSbr_dqtma(branch(il,:), V, 2);
%      [Hf16, Hf26, Hf36, Hf46, Hf56, Hf61, Hf62, Hf63, Hf64, Hf65, Hf66] = d2Abr_dxqtma2(d2Sbr_dqtma2, dSf_dVa, dSf_dVm, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dQtma, Sf, V, mu)
%      [Ht16, Ht26, Ht36, Ht46, Ht56, Ht61, Ht62, Ht63, Ht64, Ht65, Ht66] = d2Abr_dxqtma2(d2Sbr_dqtma2, dSt_dVa, dSt_dVm, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dQtma, St, V, mu)
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

[F16, F26, F36, F46, F56, F61, F62, F63, F64, F65, F66] = d2Sbr_dxqtma2(V,conj(Sbr).*mu);

d2Abr_dqtmaVa    = 2 * real( F16 + dSbr_dVa.'   * diagmu * conj(dSbr_dQtma) ); 
d2Abr_dqtmaVm    = 2 * real( F26 + dSbr_dVm.'   * diagmu * conj(dSbr_dQtma) );
d2Abr_dqtmaBeqz  = 2 * real( F36 + dSbr_dBeqz.' * diagmu * conj(dSbr_dQtma) );
d2Abr_dqtmaBeqv  = 2 * real( F46 + dSbr_dBeqv.' * diagmu * conj(dSbr_dQtma) );
d2Abr_dqtmaSh    = 2 * real( F56 + dSbr_dPfsh.' * diagmu * conj(dSbr_dQtma) );
d2Abr_dVaqtma    = 2 * real( F61 + dSbr_dQtma.' * diagmu * conj(dSbr_dVa  ) );
d2Abr_dVmqtma    = 2 * real( F62 + dSbr_dQtma.' * diagmu * conj(dSbr_dVm  ) );
d2Abr_dBeqzqtma  = 2 * real( F63 + dSbr_dQtma.' * diagmu * conj(dSbr_dBeqz) );
d2Abr_dBeqvqtma  = 2 * real( F64 + dSbr_dQtma.' * diagmu * conj(dSbr_dBeqv) );
d2Abr_dShqtma    = 2 * real( F65 + dSbr_dQtma.' * diagmu * conj(dSbr_dPfsh) );
d2Abr_dqtma2     = 2 * real( F66 + dSbr_dQtma.' * diagmu * conj(dSbr_dQtma) );
 
H16 = sparse(d2Abr_dqtmaVa);
H26 = sparse(d2Abr_dqtmaVm);
H36 = sparse(d2Abr_dqtmaBeqz);
H46 = sparse(d2Abr_dqtmaBeqv);
H56 = sparse(d2Abr_dqtmaSh);
H61 = sparse(d2Abr_dVaqtma);
H62 = sparse(d2Abr_dVmqtma);
H63 = sparse(d2Abr_dBeqzqtma);
H64 = sparse(d2Abr_dBeqvqtma);
H65 = sparse(d2Abr_dShqtma);
H66 = sparse(d2Abr_dqtma2);
