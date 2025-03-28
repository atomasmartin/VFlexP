function fun_init_intconv_DC_VSM(block, vmod, vangle, P, Q, grid_block)
    mdlWks = get_param(bdroot(grid_block), 'ModelWorkspace');
    bvarname = get_param(grid_block, 'bvarname');
    DC = evalin(mdlWks, [bvarname 'DC']);
    init_states.P = P;
    init_states.Q = Q;
    bvarname = get_param(grid_block, 'bvarname');
    block_name = get_param(block, 'Name');
    vd_out = real(vmod.*exp(vangle*1i));
    vq_out = imag(vmod.*exp(vangle*1i));
    id_out = (vd_out.*P + vq_out.*Q)./(vmod.^2);
    iq_out = (vq_out.*P - vd_out.*Q)./(vmod.^2);
    init_states.id_out = id_out;
    init_states.iq_out = iq_out;
    set_param([block '/RL out/dynamic/i_d'], 'InitialCondition', [bvarname 'op' block_name '.id_out'])
    set_param([block '/RL out/dynamic/i_q'], 'InitialCondition', [bvarname 'op' block_name '.iq_out'])
    R_out = DC.(block_name + "_dparams").Rf;
    L_out = DC.(block_name + "_dparams").Lf;
    vd_c = vd_out + R_out.*init_states.id_out - L_out.*init_states.iq_out;
    vq_c = vq_out + L_out.*init_states.id_out + R_out.*init_states.iq_out;
    init_states.frame = angle(vd_c + vq_c*1i);
    vd_in = vd_out + R_out.*init_states.id_out - L_out.*init_states.iq_out;
    vq_in = vq_out + L_out.*init_states.id_out + R_out.*init_states.iq_out;
    Kint = evalin('base', get_param(block, 'Kint'));
    Bint = evalin('base', get_param(block, 'Bint'));
    Dv = evalin('base', get_param(block, 'Dv'));
    init_states.id_in = (vd_in - vd_out + L_out.*init_states.iq_out)./Kint - (Bint.*init_states.id_out - init_states.id_out);
    init_states.iq_in = (vq_in - vq_out - L_out.*init_states.id_out)./Kint - (Bint.*init_states.iq_out - init_states.iq_out);
    set_param([block '/VCC/dynamic/PIid/I'], 'InitialCondition', [bvarname 'op' block_name '.id_in'])
    set_param([block '/VCC/dynamic/PIiq/I'], 'InitialCondition', [bvarname 'op' block_name '.iq_in'])
    init_states.IP = 0*Dv; % washout filter
    init_states.iv = sqrt(vd_in.^2 + vq_in.^2) - vmod;
    assignin(mdlWks, [bvarname 'op' block_name], init_states)
    set_param([block '/VCC/dynamic/RL/dynamic/i_d'], 'InitialCondition', [bvarname 'op' block_name '.id_in'])
    set_param([block '/VCC/dynamic/RL/dynamic/i_q'], 'InitialCondition', [bvarname 'op' block_name '.iq_in'])
    set_param([block '/frame'], 'InitialCondition', [bvarname 'op' block_name '.frame'])
    set_param([block '/PC/omega_i'], 'InitialCondition', '1')
    set_param([block '/PC/IP'], 'InitialCondition', [bvarname 'op' block_name '.IP'])
    set_param([block '/PC/LPF_P/I'], 'InitialCondition', [bvarname 'op' block_name '.P'])
    set_param([block '/PC/LPF_Q/I'], 'InitialCondition', '0')
    set_param([block '/iv'], 'InitialCondition', [bvarname 'op' block_name '.iv'])
end