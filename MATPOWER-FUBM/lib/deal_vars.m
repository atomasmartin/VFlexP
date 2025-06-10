function [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, type)
%DEAL_VARS  Deals x variables depending on the active ones.
%   [VA, VM, PG, QG, BEQZ, BEQV, SHANG, MAQT, MAVT, SHANGDP] = DEAL_VARS(X, NBEQZ, NBEQV, NPFSH, NQTMA, NVTMA, NPFDP, TYPE)
%
%   Reads the x vector of variables, and depending on which ones are active
%   it will return the adecuate value. The ammount of outputs with value 
%   will vary depending on the ammount of active variables.
%
%   Inputs:
%     X     : optimization vector
%     NBEQZ : number of elements with zero constraint control by Beq (Qf = 0)
%     NBEQV : number of elements with voltage "from" controlled by Beq    (Vf = Vf_set)
%     NPFSH : number of elements with active power "from" controlled by Theta_sh (Pf = Pf_set)
%     NQTMA : number of elements with reactive power "to" controlled by ma/tap (Qt = Qt_set)
%     NVTMA : number of elements with voltage "to" controlled by ma/tap (Vt = Vt_set)
%     NPFDP : number of elements with active power "from" controlled by Theta_sh {Pf - Pfset = Kdp*(Vmf - Vmfset)}
%     TYPE  : type of output: 
%                           "1" for nodal_balance_vars  | all vars
%                           "2" for flow_lim_vars       | all but Pg Qg
%                           "3" for zero_lim_vars       | all but Pg Qg Beqv
%                           "4" for qtma_lim_vars       | all but Pg Qg Vtma 
%                           "5" for pfsh_lim_vars       | all but Pg Qg ShAngDP
%                           "6" for pfdp_lim_vars       | all but Pg Qg ShAng
%
%   Outputs:
%     VA      : Voltage  Angle                 [rad]
%     VM      : Voltage  magnitude              [pu]
%     Pg      : Active   Generation             [pu]
%     Qg      : Reactive Generation             [pu]
%     Beqz    : Beq      for Control Qf=0       [pu]
%     Beqv    : Beq      for Control Vf=Vf_set  [pu]
%     ShAng   : Theta_sh for Control Pf=Pf_set [rad]
%     maQt    : ma/tap   for Control Qt=Qt_set  [pu]
%     maVt    : ma/tap   for Control Vt=Vt_set  [pu]
%     ShAngDP : Theta_dp for Control Pf - Pfset = Kdp*(Vmf - Vmfset) [rad]
%
%   Examples:
%       [Va, Vm, Pg, Qg, 0   , 0   , 0    , 0   , 0   , 0      ] = deal_vars(x, type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , 0    , 0   , 0   , 0      ] = deal_vars(x, nBeqz, 0    , 0    , 0    , 0    , 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, Beqv, 0    , 0   , 0   , 0      ] = deal_vars(x, nBeqz, nBeqv, 0    , 0    , 0    , 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, 0   , 0   , 0      ] = deal_vars(x, nBeqz, nBeqv, nPfsh, 0    , 0    , 0    , type); 
%       [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, 0   , 0      ] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, 0    , 0    , type); 
%       [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, 0      ] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal_vars(x, nBeqz, nBeqv, nPfsh, nQtma, nVtma, nPfdp, type);
%
%   or any other combination, ie:
%       [Va, Vm, Pg, Qg, 0   , 0   , ShAng, 0   , 0   , 0      ] = deal_vars(x, 0    , 0    , nPfsh,   0    , 0    , 0    , type); 
%       [Va, Vm, Pg, Qg, 0   , 0   , 0    , maQt, 0   , 0      ] = deal_vars(x, 0    , 0    , 0    ,   nQtma, 0    , 0    , type);
%       [Va, Vm, Pg, Qg, 0   , 0   , 0    , 0   , maVt, 0      ] = deal_vars(x, 0    , 0    , 0    ,   0    , nVtma, 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, 0   , 0   , 0      ] = deal_vars(x, nBeqz, 0    , nPfsh,   0    , 0    , 0    , type); 
%       [Va, Vm, Pg, Qg, Beqz, 0   , 0    , maQt, 0   , 0      ] = deal_vars(x, nBeqz, 0    , 0    ,   nQtma, 0    , 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, maQt, 0   , 0      ] = deal_vars(x, nBeqz, 0    , nPfsh,   nQtma, 0    , 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, maQt, 0   , 0      ] = deal_vars(x, nBeqz, 0    , nPfsh,   nQtma, 0    , 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , 0    , 0   , maVt, 0      ] = deal_vars(x, nBeqz, 0    , 0    ,   0    , nVtma, 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, 0   , maVt, 0      ] = deal_vars(x, nBeqz, 0    , nPfsh,   0    , nVtma, 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, 0   , maVt, 0      ] = deal_vars(x, nBeqz, 0    , nPfsh,   0    , nVtma, 0    , type);
%       [Va, Vm, Pg, Qg, Beqz, 0   , ShAng, 0   , maVt, ShAngDp] = deal_vars(x, nBeqz, 0    , nPfsh,   0    , nVtma, nPfdp, type);
%
%   See also OPF_POWER_BALANCE_FCN, OPF_POWER_BALANCE_HESS
                                           
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%keyboard
%% Identify amount of variables
if nargin < 2 %Only Va, Vm, Pg, Qg
    if type == 1 %nodal_balance_vars
        [Va, Vm, Pg, Qg] = deal(x{:});
        Beqz = 0;
        Beqv = 0;
        ShAng = 0;
        maQt = 0;
        maVt = 0;
    elseif type == 2 %flow_lim_vars
        [Va, Vm] = deal(x{:});
        Pg = 0;
        Qg = 0;
        Beqz = 0;
        Beqv = 0; 
        ShAng = 0;
        maQt = 0;
        maVt = 0;
    else
        error('deal_vars: type can only be values from 1 to 2 for less than 2 inputs')
    end
else %FUBM variables included
    if nBeqz
        if nBeqv 
            if nPfsh
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal(x{:});
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, maVt] = deal(x{:});
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt, ShAngDp] = deal(x{:});
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maQt] = deal(x{:});
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maVt, ShAngDp] = deal(x{:});
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, maVt] = deal(x{:});
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng, ShAngDp] = deal(x{:});
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAng] = deal(x{:});
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            else %nPfsh = 0
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:});
                                ShAng = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maQt, maVt] = deal(x{:});
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maQt, ShAngDp] = deal(x{:});
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maQt] = deal(x{:});
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maVt, ShAngDp] = deal(x{:});
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, maVt] = deal(x{:});
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv, ShAngDp] = deal(x{:});
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, Beqv] = deal(x{:});
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            end
        else % nBeqv = 0
            if nPfsh
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maQt, maVt] = deal(x{:});
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maQt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maQt] = deal(x{:});
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maVt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, maVt] = deal(x{:});
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAng] = deal(x{:});
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            else %nPfsh = 0
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maQt, maVt] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maQt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maQt] = deal(x{:});
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maVt, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, maVt] = deal(x{:});
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz, ShAngDp] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqz] = deal(x{:});
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqz] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            end
        end
    else %nBeqz = 0
        if nBeqv 
            if nPfsh
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maQt, maVt] = deal(x{:});
                                Beqz = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maQt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maQt] = deal(x{:});
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, maVt] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Beqz = 0;
                                Pg = 0;
                                Qg = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAng] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            else %nPfsh = 0
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maQt, maVt] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maQt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maQt] = deal(x{:});
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, maVt] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, Beqv] = deal(x{:});
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, Beqv] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            end
        else % nBeqv = 0
            if nPfsh
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maQt, maVt] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maQt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maQt] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Beqz = 0;
                                Pg = 0;
                                Qg = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, maVt] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAng] = deal(x{:});
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, ShAng] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            else %nPfsh = 0
                if nQtma
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maQt, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maQt, maVt] = deal(x{:});
                                Beqv = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maQt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, maQt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maQt] = deal(x{:});
                                Beqv = 0;
                                maVt = 0;
                                Beqz = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maQt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maVt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                else %nQtma = 0
                    if nVtma
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maVt, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maVt = 0;
                                maQt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, maVt] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                maQt = 0;
                                ShAng = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    else %nVtma = 0;
                        if nPfdp 
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg, ShAngDp] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm, ShAngDp] = deal(x{:}); 
                                Beqz = 0;
                                Pg = 0;
                                Qg = 0;
                                ShAng = 0;
                                maQt = 0;
                                Beqv = 0;
                                maVt = 0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm, maVt, ShAngDp] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        else %nPfdp = 0;
                            if type == 1 %nodal_balance_vars | all vars
                                [Va, Vm, Pg, Qg] = deal(x{:});
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 2 %flow_lim_vars  | all but Pg Qg
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 3 %zero_lim_vars  | all but Pg Qg Beqv
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                Beqv = 0;
                                ShAngDp=0;
                            elseif type == 4 %qtma_lim_vars  | all but Pg Qg Vtma 
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            elseif type == 5 %pfsh_lim_vars  | all but Pg Qg ShAngDP
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp = 0;
                            elseif type == 6 %pfdp_lim_vars  | all but Pg Qg ShAng
                                [Va, Vm] = deal(x{:}); 
                                Pg = 0;
                                Qg = 0;
                                Beqz = 0;
                                Beqv = 0;
                                ShAng = 0;
                                maQt = 0;
                                maVt = 0;
                                ShAngDp=0;
                            else
                                error('deal_vars: type can only be values from 1 to 6') 
                            end 
                        end
                    end
                end
            end
        end
    end
end