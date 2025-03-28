function fun_init_GFr(block, vmod, vangle, P, Q, grid_block)
    mdlWks = get_param(bdroot(grid_block), 'ModelWorkspace');
    set_param([block '/PI_synchro_f/I'], 'InitialCondition', '0')
    set_param([block '/PI_synchro_v/I'], 'InitialCondition', '0')
    set_param([block '/LPF_v/I'], 'InitialCondition', '0')
    set_param([block '/LPF_omega/I'], 'InitialCondition', '1')
    vd_out = real(vmod.*exp(vangle*1i));
    vq_out = imag(vmod.*exp(vangle*1i));
    init_states.Pd = P;
    init_states.Qd = Q;
    init_states.id_out = (vd_out.*P + vq_out.*Q)./(vmod.^2);
    init_states.iq_out = (vq_out.*P - vd_out.*Q)./(vmod.^2);
    R_out = evalin('base', get_param(block, 'R_out'));
    L_out = evalin('base', get_param(block, 'L_out'));
    Rcf = evalin('base', get_param(block, 'Rcf'));
    Cfi = evalin('base', get_param(block, 'Cfi'));
    Rf = evalin('base', get_param(block, 'Rf'));
    Lf = evalin('base', get_param(block, 'Lf'));
    Kint = evalin('base', get_param(block, 'Kint'));
    Bint = evalin('base', get_param(block, 'Bint'));
    Tint = evalin('base', get_param(block, 'Tint'));
    Kext = evalin('base', get_param(block, 'Kext'));
    Bext = evalin('base', get_param(block, 'Bext'));
    Text = evalin('base', get_param(block, 'Text'));
    Fi = evalin('base', get_param(block, 'Fi'));
    vd_c = vd_out + R_out.*init_states.id_out - L_out.*init_states.iq_out;
    vq_c = vq_out + L_out.*init_states.id_out + R_out.*init_states.iq_out;
    init_states.frame = angle(vd_c + vq_c*1i);
    init_states.vd_c = real(exp(1i*-init_states.frame).*(vd_c + 1i*vq_c));
    init_states.vq_c = imag(exp(1i*-init_states.frame).*(vd_c + 1i*vq_c));
    vd_out_own = real(exp(1i*-init_states.frame).*(vd_out + 1i*vq_out));
    vq_out_own = imag(exp(1i*-init_states.frame).*(vd_out + 1i*vq_out));
    init_states.iV = init_states.vd_c - vmod;
    init_states.v = init_states.vd_c;
    init_states.IntConnect = 5*ones(1, length(init_states.vd_c));
    id_out_own = real(exp(1i*-init_states.frame).*(init_states.id_out + 1i*init_states.iq_out));
    iq_out_own = imag(exp(1i*-init_states.frame).*(init_states.id_out + 1i*init_states.iq_out));
    init_states.id_l = id_out_own + init_states.vd_c./Rcf - Cfi.*init_states.vq_c;
    init_states.iq_l = iq_out_own + Cfi.*init_states.vd_c + init_states.vq_c./Rcf;
    bvarname = get_param(grid_block, 'bvarname');
    block_name = get_param(block, 'Name');
    init_states.PIid = Rf.*init_states.id_l/Kint - (Bint.*init_states.id_l - init_states.id_l);
    init_states.PIiq = Rf.*init_states.iq_l/Kint - (Bint.*init_states.iq_l - init_states.iq_l);
    init_states.PIvd = (init_states.id_l - Fi.*id_out_own + Cfi.*init_states.vq_c)/Kext - (Bext.*init_states.vd_c - init_states.vd_c);
    init_states.PIvq = (init_states.iq_l - Fi.*iq_out_own - Cfi.*init_states.vd_c)/Kext - (Bext.*init_states.vq_c - init_states.vq_c);
    assignin(mdlWks, [bvarname 'op' block_name], init_states)
    set_param([block '/LPF_v/I'], 'InitialCondition', [bvarname 'op' block_name '.v'])
    set_param([block '/VCC/dynamic/RC/dynamic/v_d'], 'InitialCondition', [bvarname 'op' block_name '.vd_c'])
    set_param([block '/VCC/dynamic/RC/dynamic/v_q'], 'InitialCondition', [bvarname 'op' block_name '.vq_c'])
    set_param([block '/VCC/dynamic/RL/dynamic/i_d'], 'InitialCondition', [bvarname 'op' block_name '.id_l'])
    set_param([block '/VCC/dynamic/RL/dynamic/i_q'], 'InitialCondition', [bvarname 'op' block_name '.iq_l'])
    set_param([block '/VCC/dynamic/PIvd/I'], 'InitialCondition', [bvarname 'op' block_name '.PIvd'])
    set_param([block '/VCC/dynamic/PIvq/I'], 'InitialCondition', [bvarname 'op' block_name '.PIvq'])
    set_param([block '/VCC/dynamic/PIid/I'], 'InitialCondition', [bvarname 'op' block_name '.PIid'])
    set_param([block '/VCC/dynamic/PIiq/I'], 'InitialCondition', [bvarname 'op' block_name '.PIiq'])
    set_param([block '/RL out/dynamic/i_d'], 'InitialCondition', [bvarname 'op' block_name '.id_out'])
    set_param([block '/RL out/dynamic/i_q'], 'InitialCondition', [bvarname 'op' block_name '.iq_out'])
    set_param([block '/frame'], 'InitialCondition', [bvarname 'op' block_name '.frame'])
    set_param([block '/PC/LPF_P/I'], 'InitialCondition', '0')
    set_param([block '/PC/LPF_Q/I'], 'InitialCondition', '0')
    set_param([block '/iV'], 'InitialCondition', [bvarname 'op' block_name '.iV'])
    set_param([block '/Integrator'], 'InitialCondition', [bvarname 'op' block_name '.IntConnect'])
end