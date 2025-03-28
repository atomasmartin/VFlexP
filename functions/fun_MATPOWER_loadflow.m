function mpc_result = fun_MATPOWER_loadflow(gcb, writebranch)
    try
        if nargin == 1
            writebranch = 0;
        end
        mpc = fun_VFlexP_to_MATPOWER(gcb);
        if isempty(mpc.branch)
            mpoption0 = mpoption;
            mpoption0.out.all = 0;
            mpc_result = runpf(mpc, mpoption0);
            fun_MATPOWER_to_VFlexP(mpc_result, gcb, 0)
            bvarname = get_param(gcb, 'bvarname');
            intconv_PF.include = false;
            mdlWks = get_param(bdroot(gcb), 'ModelWorkspace');
            assignin(mdlWks, [bvarname 'intconv_PF'], intconv_PF)
            fun_initialise_system_from_tables(gcb)
            fun_apply_changes_block(gcb)
            fun_init_units_callback(gcb)
        else
            mpc_result = runpf(mpc);
            mpc_result.branch_complete = mpc_result.branch;
            mpc_result.branch = mpc_result.branch(:,1:13);
            intconv_PF.include = false;
            mdlWks = get_param(bdroot(gcb), 'ModelWorkspace');
            if mpc_result.success == 1
                if size(mpc_result.branch_complete, 2) == 37 % FUBM formulation
                    idx_intconv = find(mpc_result.branch_complete(:,26) ~= 0); % interface converters
                    Nintconv = length(idx_intconv);
                    intconv_PF.vmod = zeros(1, Nintconv);
                    intconv_PF.vangle = intconv_PF.vmod;
                    intconv_PF.P = intconv_PF.vmod;
                    intconv_PF.Q = intconv_PF.vmod;
                    %from_intconv = mpc_result.branch_complete(idx_intconv,1);
                    to_intconv = mpc_result.branch_complete(idx_intconv,2);
                    bus_i = mpc_result.bus(:,1);
                    bus_vmod = mpc_result.bus(:,8);
                    bus_vangle = mpc_result.bus(:,9);
                    intconv_PF.from = mpc_result.branch(:, 1)';
                    intconv_PF.to = mpc_result.branch(:, 2)';
                    for i=1:Nintconv                   
                        intconv_PF.vmod(i) = bus_vmod(find(bus_i == to_intconv(i), 1));
                        intconv_PF.vangle(i) = bus_vangle(find(bus_i == to_intconv(i), 1));
                        intconv_PF.P(i) = mpc_result.branch_complete(idx_intconv(i), 16)/mpc.baseMVA;
                        intconv_PF.Q(i) = mpc_result.branch_complete(idx_intconv(i), 17)/mpc.baseMVA;
                    end
                    intconv_PF.include = true;
                    if writebranch == 0
                        fun_MATPOWER_to_VFlexP(mpc_result, gcb, 4)
                    else
                        fun_MATPOWER_to_VFlexP(mpc_result, gcb, 2)
                    end
                else
                    if writebranch == 0
                        fun_MATPOWER_to_VFlexP(mpc_result, gcb, 3)
                    else
                        fun_MATPOWER_to_VFlexP(mpc_result, gcb, 0)
                    end
                end
                bvarname = get_param(gcb, 'bvarname');
                assignin(mdlWks, [bvarname 'intconv_PF'], intconv_PF)
                fun_initialise_system_from_tables(gcb)
                fun_apply_changes_block(gcb)
                % fun_init_units_callback(gcb) %% do not initialise every
                % loadflow
            else
                errordlg("MATPOWER PF did not succeed.")
            end
        end
    catch ME
        str_name = string({ME.stack.name}');
        str_line = string({ME.stack.line}');
        errortext = [ME.message; strcat("Error in file: ", str_name, ": line ", str_line)];
        errordlg(char(errortext));
    end
end
    