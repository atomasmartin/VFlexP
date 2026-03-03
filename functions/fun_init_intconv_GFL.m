function fun_init_intconv_DC_GFl(block, vmod, vangle, P, Q, grid_block)
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
    init_states.id_out = real((id_out+1i*iq_out).*exp(-vangle*1i));
    init_states.iq_out = imag((id_out+1i*iq_out).*exp(-vangle*1i));
    R_out = DC.(block_name + "_dparams").Rf;
    L_out = DC.(block_name + "_dparams").Lf;
    vd_in = vd_out + R_out.*init_states.id_out - L_out.*init_states.iq_out;
    vq_in = vq_out + L_out.*init_states.id_out + R_out.*init_states.iq_out;
    Kint = evalin('base', get_param(block, 'Kint'));
    Bint = evalin('base', get_param(block, 'Bint'));
    init_states.PId = (vd_in - vd_out + L_out.*init_states.iq_out)./Kint - (Bint.*init_states.id_out - init_states.id_out);
    init_states.PIq = (vq_in - vq_out - L_out.*init_states.id_out)./Kint - (Bint.*init_states.iq_out - init_states.iq_out);
    set_param([block '/CC/dynamic/PId/I'], 'InitialCondition', [bvarname 'op' block_name '.PId'])
    set_param([block '/CC/dynamic/PIq/I'], 'InitialCondition', [bvarname 'op' block_name '.PIq'])
    init_states.frame = vangle;
    assignin(mdlWks, [bvarname 'op' block_name], init_states)
    set_param([block '/CC/dynamic/RL/dynamic/i_d'], 'InitialCondition', [bvarname 'op' block_name '.id_out'])
    set_param([block '/CC/dynamic/RL/dynamic/i_q'], 'InitialCondition', [bvarname 'op' block_name '.iq_out'])
    set_param([block '/PLL/dynamic/frame'], 'InitialCondition', [bvarname 'op' block_name '.frame'])
    set_param([block '/PLL/dynamic/PLL/I'], 'InitialCondition', '1')
    set_param([block '/P'], 'InitialCondition', [bvarname 'op' block_name '.P'])
    set_param([block '/I'], 'InitialCondition', [bvarname 'op' block_name '.P'])
    set_param([block '/Q'], 'InitialCondition', [bvarname 'op' block_name '.Q'])
end