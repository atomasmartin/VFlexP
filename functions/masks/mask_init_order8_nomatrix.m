function Matrix_flux_to_i = mask_init_order8_nomatrix(Lssd_pu, Lssq_pu, Lmd_pu, Lmq_pu, Lfd_pu, Lkd1_pu, Lkd2_pu, Lkq1_pu, Lkq2_pu, Lkq3_pu)
%MASK_INIT_ORDER8_NOMATRIX Compute flux-to-current matrix for 8th-order synchronous machine.

n = numel(Lssd_pu);
Lssd = Lssd_pu;
Lssq = Lssq_pu;
Lmd  = Lmd_pu;
Lmq  = Lmq_pu;
Lfd  = Lfd_pu;
Lkd1 = Lkd1_pu;
Lkd2 = Lkd2_pu;
Lkq1 = Lkq1_pu;
Lkq2 = Lkq2_pu;
Lkq3 = Lkq3_pu;

Matrix_flux_to_i = zeros(8, 8, max(1, n));
for i = 1:n
    dd = 2*Lmd(i)^4 - Lkd1(i)*Lmd(i)^3 - 2*Lkd2(i)*Lmd(i)^3 - Lmd(i)^3*Lssd(i) - 2*Lfd(i)*Lmd(i)^3 ...
        + Lfd(i)*Lkd1(i)*Lmd(i)^2 + Lfd(i)*Lkd2(i)*Lmd(i)^2 + Lkd1(i)*Lkd2(i)*Lmd(i)^2 ...
        + Lfd(i)*Lmd(i)^2*Lssd(i) + Lkd2(i)*Lmd(i)^2*Lssd(i) - Lfd(i)*Lkd1(i)*Lkd2(i)*Lssd(i);
    dq = 3*Lmq(i)^4 - 2*Lkq2(i)*Lmq(i)^3 - 2*Lkq3(i)*Lmq(i)^3 - 2*Lmq(i)^3*Lssq(i) - 2*Lkq1(i)*Lmq(i)^3 ...
        + Lkq1(i)*Lkq2(i)*Lmq(i)^2 + Lkq1(i)*Lkq3(i)*Lmq(i)^2 + Lkq2(i)*Lkq3(i)*Lmq(i)^2 ...
        + Lkq1(i)*Lmq(i)^2*Lssq(i) + Lkq2(i)*Lmq(i)^2*Lssq(i) + Lkq3(i)*Lmq(i)^2*Lssq(i) ...
        - Lkq1(i)*Lkq2(i)*Lkq3(i)*Lssq(i);

    M = zeros(8, 8);
    M(1,1) = -(Lfd(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 - Lmd(i)^3 - Lfd(i)*Lkd1(i)*Lkd2(i)) / dd;
    M(1,3) =  (Lkd1(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 - Lmd(i)^3 - Lkd1(i)*Lkd2(i)*Lmd(i)) / dd;
    M(1,4) =  (Lfd(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 - Lmd(i)^3 - Lfd(i)*Lkd2(i)*Lmd(i)) / dd;
    M(1,5) =  (Lfd(i)*Lmd(i)^2 - Lfd(i)*Lkd1(i)*Lmd(i)) / dd;

    M(2,2) = -(Lkq1(i)*Lmq(i)^2 + Lkq2(i)*Lmq(i)^2 + Lkq3(i)*Lmq(i)^2 - 2*Lmq(i)^3 - Lkq1(i)*Lkq2(i)*Lkq3(i)) / dq;
    M(2,6) = -(3*Lmq(i)^3 - 2*Lkq2(i)*Lmq(i)^2 - 2*Lkq3(i)*Lmq(i)^2 - 2*Lkq1(i)*Lmq(i)^2 ...
               + Lkq1(i)*Lkq2(i)*Lmq(i) + Lkq1(i)*Lkq3(i)*Lmq(i) + Lkq2(i)*Lkq3(i)*Lmq(i)) / dq;

    M(3,1) = -(Lkd2(i)*Lmd(i)^2 - Lkd1(i)*Lkd2(i)*Lmd(i)) / dd;
    M(3,3) =  (Lkd1(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - 2*Lmd(i)^3 - Lkd1(i)*Lkd2(i)*Lssd(i)) / dd;
    M(3,4) = -(Lkd2(i)*Lmd(i)^2 - Lkd2(i)*Lssd(i)*Lmd(i)) / dd;
    M(3,5) = -(Lkd1(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - 2*Lmd(i)^3) / dd;

    M(4,1) = -(Lfd(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 - Lmd(i)^3 - Lfd(i)*Lkd2(i)*Lmd(i)) / dd;
    M(4,3) = -(Lkd2(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - Lmd(i)^3 - Lkd2(i)*Lmd(i)*Lssd(i)) / dd;
    M(4,4) =  (Lfd(i)*Lmd(i)^2 + Lkd2(i)*Lmd(i)^2 - Lmd(i)^3 - Lfd(i)*Lkd2(i)*Lssd(i)) / dd;
    M(4,5) = -(Lfd(i)*Lmd(i)^2 - Lfd(i)*Lssd(i)*Lmd(i)) / dd;

    M(5,1) = -(Lfd(i)*Lmd(i)^2 + Lkd1(i)*Lmd(i)^2 - Lmd(i)^3 - Lfd(i)*Lkd1(i)*Lmd(i)) / dd;
    M(5,3) = -(Lkd1(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - Lmd(i)^3 - Lkd1(i)*Lmd(i)*Lssd(i)) / dd;
    M(5,4) = -(Lfd(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - Lmd(i)^3 - Lfd(i)*Lmd(i)*Lssd(i)) / dd;
    M(5,5) =  (Lfd(i)*Lmd(i)^2 + Lkd1(i)*Lmd(i)^2 + Lmd(i)^2*Lssd(i) - 2*Lmd(i)^3 - Lfd(i)*Lkd1(i)*Lssd(i)) / dd;

    M(6,2) = -(Lkq2(i)*Lmq(i)^2 + Lkq3(i)*Lmq(i)^2 - Lmq(i)^3 - Lkq2(i)*Lkq3(i)*Lmq(i)) / dq;
    M(6,6) = -(Lmq(i)^2*Lssq(i) + Lkq2(i)*Lkq3(i)*Lssq(i) - Lkq2(i)*Lmq(i)*Lssq(i) - Lkq3(i)*Lmq(i)*Lssq(i)) / dq;

    M(7,2) =  (Lmq(i)*(Lkq1(i) - Lmq(i))*(Lkq3(i) - Lmq(i))) / dq;
    M(7,6) = -(Lssq(i)*(Lkq1(i) - Lmq(i))*(Lkq3(i) - Lmq(i))) / dq;

    M(8,2) =  (Lmq(i)*(Lkq1(i) - Lmq(i))*(Lkq2(i) - Lmq(i))) / dq;
    M(8,6) = -(Lssq(i)*(Lkq1(i) - Lmq(i))*(Lkq2(i) - Lmq(i))) / dq;

    Matrix_flux_to_i(:,:,i) = M;
end
end
