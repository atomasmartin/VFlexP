function [H15, H25, H35, H45, H51, H52, H53, H54, H55] = d2Abr_dxsh2(d2Sbr_dxsh2, dSbr_dVa, dSbr_dVm, dSbr_dBeqz, dSbr_dBeqv, dSbr_dPfsh, Sbr, V, mu)
%D2ABR_DXSH2  Computes 2nd derivatives of |branch flow|^2 w.r.t. ShVa, ShVm, ShBeqz, ShBeqv, VaSh, VmSh, BeqzSh, BeqvSh, ShSh.
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H15, H25, H35, H45, H51, H52, H53, H54, H55] = D2ABR_DXSH2(D2SBR_DXSH2, DSBR_DVA, DSBR_DVM, DSBR_DBEQZ, DSBR_DBEQV, DSBR_DPFSH, SBR, V, MU)
%
%   Returns 9 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H15 = HShVa = d2Abr_dShVa     = d/dsh   ( (dAbr/dVa  )'*mu )
%   H25 = HShVm = d2Abr_dShVm     = d/dsh   ( (dAbr/dVm  )'*mu )
%   H35 = HShBz = d2Abr_dShBeqz   = d/dsh   ( (dAbr/dBeqz)'*mu )
%   H45 = HShBv = d2Abr_dShBeqv   = d/dsh   ( (dAbr/dBeqv)'*mu )
%   H51 = HVaSh = d2Abr_dVaSh     = d/dVa   ( (dAbr/dsh  )'*mu )
%   H52 = HVmSh = d2Abr_dVmSh     = d/dVm   ( (dAbr/dsh  )'*mu )
%   H53 = HBzSh = d2Abr_dBeqzSh   = d/dBeqz ( (dAbr/dsh  )'*mu )
%   H54 = HBvSh = d2Abr_dBeqvSh   = d/dBeqv ( (dAbr/dsh  )'*mu )
%   H55 = HShSh = d2Abr_dShSh     = d/dsh   ( (dAbr/dsh  )'*mu )
%   
%   H15 = HShVa = d2Abr_dShVa     =  (  2*real(d2Sbr_dShVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dSh   ))'*mu 
%   H25 = HShVm = d2Abr_dShVm     =  (  2*real(d2Sbr_dShVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dSh   ))'*mu
%   H35 = HShBz = d2Abr_dShBeqz   =  (  2*real(d2Sbr_dShBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dSh   ))'*mu
%   H45 = HShBv = d2Abr_dShBeqv   =  (  2*real(d2Sbr_dShBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dSh   ))'*mu
%   H51 = HVaSh = d2Abr_dVaSh     =  (  2*real(d2Sbr_dVaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dVa   ))'*mu
%   H52 = HVmSh = d2Abr_dVmSh     =  (  2*real(d2Sbr_dVmSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dVm   ))'*mu
%   H53 = HBzSh = d2Abr_dBeqzSh   =  (  2*real(d2Sbr_dBeqzSh)*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dBeqz ))'*mu
%   H54 = HBvSh = d2Abr_dBeqvSh   =  (  2*real(d2Sbr_dBeqvSh)*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dBeqv ))'*mu
%   H55 = HShSh = d2Abr_dShSh     =  (  2*real(d2Sbr_dShSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dSh   ))'*mu 
%
%   Example:
%
%      [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch(il,:));
%      d2Sf_dsh2 = @(V, mu)d2Sf_dxsh2(branch(il,:), V, mu, mpopt.opf.v_cartesian);
%      [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%      [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch(il,:), V, 1);
%      [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch(il,:), V, 2);
%      [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch(il,:), V, 1);
%      [Hf15, Hf25, Hf35, Hf45, Hf51, Hf52, Hf53, Hf54, Hf55] = d2Abr_dxsh2(d2Sbr_dsh2, dSf_dVa, dSf_dVm, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, Sf, V, mu)
%      [Ht15, Ht25, Ht35, Ht45, Ht51, Ht52, Ht53, Ht54, Ht55] = d2Abr_dxsh2(d2Sbr_dsh2, dSt_dVa, dSt_dVm, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, St, V, mu)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV, D2ABR_DV2, D2ABR_DXBEQV2.
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

[F15, F25, F35, F45, F51, F52, F53, F54, F55] = d2Sbr_dxsh2(V,conj(Sbr).*mu);

d2Abr_dShVa    = 2 * real( F15 + dSbr_dVa.'   * diagmu * conj(dSbr_dPfsh) ); 
d2Abr_dShVm    = 2 * real( F25 + dSbr_dVm.'   * diagmu * conj(dSbr_dPfsh) );
d2Abr_dShBeqz  = 2 * real( F35 + dSbr_dBeqz.' * diagmu * conj(dSbr_dPfsh) );
d2Abr_dShBeqv  = 2 * real( F45 + dSbr_dBeqv.' * diagmu * conj(dSbr_dPfsh) );
d2Abr_dVaSh    = 2 * real( F51 + dSbr_dPfsh.' * diagmu * conj(dSbr_dVa  ) );
d2Abr_dVmSh    = 2 * real( F52 + dSbr_dPfsh.' * diagmu * conj(dSbr_dVm  ) );
d2Abr_dBeqzSh  = 2 * real( F53 + dSbr_dPfsh.' * diagmu * conj(dSbr_dBeqz) );
d2Abr_dBeqvSh  = 2 * real( F54 + dSbr_dPfsh.' * diagmu * conj(dSbr_dBeqv) );
d2Abr_dSh2     = 2 * real( F55 + dSbr_dPfsh.' * diagmu * conj(dSbr_dPfsh) );
 
H15 = sparse(d2Abr_dShVa);
H25 = sparse(d2Abr_dShVm);
H35 = sparse(d2Abr_dShBeqz);
H45 = sparse(d2Abr_dShBeqv);
H51 = sparse(d2Abr_dVaSh);
H52 = sparse(d2Abr_dVmSh);
H53 = sparse(d2Abr_dBeqzSh);
H54 = sparse(d2Abr_dBeqvSh);
H55 = sparse(d2Abr_dSh2);
