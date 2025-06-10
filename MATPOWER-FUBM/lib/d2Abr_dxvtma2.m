function [H17, H27, H37, H47, H57, H67, H71, H72, H73, H74, H75, H76, H77] = d2Abr_dxvtma2(d2Sbr_dxvtma2, dSbr_dVa, dSbr_dVm, dSbr_dBeqz, dSbr_dBeqv, dSbr_dPfsh, dSbr_dqtma, dSbr_dvtma, Sbr, V, mu)
%D2ABR_DXVTMA2  Computes 2nd derivatives of |branch flow|^2 w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma, vtmavtma.
%
%   The derivatives will be taken with respect to polar or cartesian
%   coordinates (So far only Polar has been coded) 
%   Flows could be complex current or complex or real power. 
%  (So far only complex power has been coded). Notation below is based on 
%   complex power.
%
%   [H17, H27, H37, H47, H57, H67, H71, H72, H73, H74, H75, H76, H77] = D2ABR_DXVTMA2(D2SBR_DXVTMA2, DSBR_DVA, DSBR_DVM, DSBR_DBEQZ, DSBR_DBEQV, DSBR_DPFSH, DSBR_DQTMA, DSBR_DVTMA, SBR, V, MU)
%
%   Returns 13 matrices containing the 2nd partial derivatives of Af
%   where:
%   Abr = abs(Sbr).^2
%      = Sbr.*conj(Sbr)
%     
%   dAbr/dx = 2*real( diag(conj(Sbr))*dSbr_dx)
%
%   thus,
%
%   H17 = HvtmaVa   = d2Abr_dvtmaVa     = d/dvtma ( (dAbr/dVa    )'*mu )
%   H27 = HvtmaVm   = d2Abr_dvtmaVm     = d/dvtma ( (dAbr/dVm    )'*mu )
%   H37 = HvtmaBz   = d2Abr_dvtmaBeqz   = d/dvtma ( (dAbr/dBeqz  )'*mu )
%   H47 = HvtmaBv   = d2Abr_dvtmaBeqv   = d/dvtma ( (dAbr/dBeqv  )'*mu )
%   H57 = HvtmaSh   = d2Abr_dvtmaSh     = d/dvtma ( (dAbr/dSh    )'*mu )
%   H57 = HvtmaSh   = d2Abr_dvtmaqtma   = d/dvtma ( (dAbr/dqtma  )'*mu )
%   H71 = HVavtma   = d2Abr_dVavtma     = d/dVa   ( (dAbr/dvtma  )'*mu )
%   H72 = HVmvtma   = d2Abr_dVmvtma     = d/dVm   ( (dAbr/dvtma  )'*mu )
%   H73 = HBzvtma   = d2Abr_dBeqzvtma   = d/dBeqz ( (dAbr/dvtma  )'*mu )
%   H74 = HBvvtma   = d2Abr_dBeqvvtma   = d/dBeqv ( (dAbr/dvtma  )'*mu )
%   H75 = HShvtma   = d2Abr_dShvtma     = d/dsh   ( (dAbr/dvtma  )'*mu )
%   H76 = Hqtmavtma = d2Abr_dqtmavtma   = d/dqtma ( (dAbr/dvtma  )'*mu )
%   H77 = Hvtmavtma = d2Abr_dvtmavtma   = d/dvtma ( (dAbr/dvtma  )'*mu )
%   
%   H17 = HvtmaVa   = d2Abr_dvtmaVa     =  (  2*real(d2Sbr_dvtmaVa  )*conj(Sbr)  +  dSbr_dVa   *conj(dSbr_dvtma   ))'*mu 
%   H27 = HvtmaVm   = d2Abr_dvtmaVm     =  (  2*real(d2Sbr_dvtmaVm  )*conj(Sbr)  +  dSbr_dVm   *conj(dSbr_dvtma   ))'*mu
%   H37 = HvtmaBz   = d2Abr_dvtmaBeqz   =  (  2*real(d2Sbr_dvtmaBeqz)*conj(Sbr)  +  dSbr_dBeqz *conj(dSbr_dvtma   ))'*mu
%   H47 = HvtmaBv   = d2Abr_dvtmaBeqv   =  (  2*real(d2Sbr_dvtmaBeqv)*conj(Sbr)  +  dSbr_dBeqv *conj(dSbr_dvtma   ))'*mu
%   H57 = HvtmaSh   = d2Abr_dvtmaSh     =  (  2*real(d2Sbr_dvtmaSh  )*conj(Sbr)  +  dSbr_dSh   *conj(dSbr_dvtma   ))'*mu
%   H67 = Hvtmaqtma = d2Abr_dvtmaqtma   =  (  2*real(d2Sbr_dvtmaqtma)*conj(Sbr)  +  dSbr_dqtma *conj(dSbr_dvtma   ))'*mu
%   H71 = HVavtma   = d2Abr_dVavtma     =  (  2*real(d2Sbr_dVavtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dVa     ))'*mu
%   H72 = HVmvtma   = d2Abr_dVmvtma     =  (  2*real(d2Sbr_dVmvtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dVm     ))'*mu
%   H73 = HBzvtma   = d2Abr_dBeqzvtma   =  (  2*real(d2Sbr_dBeqzvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dBeqz   ))'*mu
%   H74 = HBvvtma   = d2Abr_dBeqvvtma   =  (  2*real(d2Sbr_dBeqvvtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dBeqv   ))'*mu
%   H75 = HShvtma   = d2Abr_dShvtma     =  (  2*real(d2Sbr_dShvtma  )*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dSh     ))'*mu 
%   H76 = Hqtmavtma = d2Abr_dqtmavtma   =  (  2*real(d2Sbr_dqtmavtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dqtma   ))'*mu 
%   H77 = Hvtmavtma = d2Abr_dvtmavtma   =  (  2*real(d2Sbr_dvtmavtma)*conj(Sbr)  +  dSbr_dvtma *conj(dSbr_dvtma   ))'*mu 
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
%      [dSf_dVtma, dSt_dVtma] = dSbr_dvtma(branch(il,:), V, 4);
%      [Hf17, Hf27, Hf37, Hf47, Hf57, Hf67, Hf71, Hf72, Hf73, Hf74, Hf75, Hf76, Hf77] = d2Abr_dxvtma2(d2Sf_dxvtma2, dSf_dVa, dSf_dVm, dSf_dBeqz, dSf_dBeqv, dSf_dPfsh, dSf_dqtma, dSf_dvtma, Sf, V, mu)
%      [Ht17, Ht27, Ht37, Ht47, Ht57, Ht67, Ht71, Ht72, Ht73, Ht74, Ht75, Ht76, Ht77] = d2Abr_dxvtma2(d2St_dxvtma2, dSt_dVa, dSt_dVm, dSt_dBeqz, dSt_dBeqv, dSt_dPfsh, dSt_dqtma, dSt_dvtma, St, V, mu)
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

[F17, F27, F37, F47, F57, F67, F71, F72, F73, F74, F75, F76, F77] = d2Sbr_dxvtma2(V,conj(Sbr).*mu);

d2Abr_dvtmaVa    = 2 * real( F17 + dSbr_dVa.'   * diagmu * conj(dSbr_dvtma) ); 
d2Abr_dvtmaVm    = 2 * real( F27 + dSbr_dVm.'   * diagmu * conj(dSbr_dvtma) );
d2Abr_dvtmaBeqz  = 2 * real( F37 + dSbr_dBeqz.' * diagmu * conj(dSbr_dvtma) );
d2Abr_dvtmaBeqv  = 2 * real( F47 + dSbr_dBeqv.' * diagmu * conj(dSbr_dvtma) );
d2Abr_dvtmaSh    = 2 * real( F57 + dSbr_dPfsh.' * diagmu * conj(dSbr_dvtma) );
d2Abr_dvtmaqtma  = 2 * real( F67 + dSbr_dqtma.' * diagmu * conj(dSbr_dvtma) );
d2Abr_dVavtma    = 2 * real( F71 + dSbr_dvtma.' * diagmu * conj(dSbr_dVa  ) );
d2Abr_dVmvtma    = 2 * real( F72 + dSbr_dvtma.' * diagmu * conj(dSbr_dVm  ) );
d2Abr_dBeqzvtma  = 2 * real( F73 + dSbr_dvtma.' * diagmu * conj(dSbr_dBeqz) );
d2Abr_dBeqvvtma  = 2 * real( F74 + dSbr_dvtma.' * diagmu * conj(dSbr_dBeqv) );
d2Abr_dShvtma    = 2 * real( F75 + dSbr_dvtma.' * diagmu * conj(dSbr_dPfsh) );
d2Abr_dqtmavtma  = 2 * real( F76 + dSbr_dvtma.' * diagmu * conj(dSbr_dqtma) );
d2Abr_dvtma2     = 2 * real( F77 + dSbr_dvtma.' * diagmu * conj(dSbr_dvtma) );
 
H17 = sparse(d2Abr_dvtmaVa);
H27 = sparse(d2Abr_dvtmaVm);
H37 = sparse(d2Abr_dvtmaBeqz);
H47 = sparse(d2Abr_dvtmaBeqv);
H57 = sparse(d2Abr_dvtmaSh);
H67 = sparse(d2Abr_dvtmaqtma);
H71 = sparse(d2Abr_dVavtma);
H72 = sparse(d2Abr_dVmvtma);
H73 = sparse(d2Abr_dBeqzvtma);
H74 = sparse(d2Abr_dBeqvvtma);
H75 = sparse(d2Abr_dShvtma);
H76 = sparse(d2Abr_dqtmavtma);
H77 = sparse(d2Abr_dvtma2);
