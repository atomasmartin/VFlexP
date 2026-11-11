function [ref, sentSC, sentAll, following_o, sending_o, delivered] = fun_MAS_tc_chart(elec, vmod, MAS_cfg_vec, MAS_Pnom, MAS_vs, MAS_isdc, MAS_buses) %#codegen
% fun_MAS_tc_chart - Vectorised multi-agent tertiary control (variant "MAS").
%
% Source file for the MATLAB Function block
%   <model>/Grid/Tertiary Control/MAS/MAS_agents
% (the block Script is a verbatim copy of this file; edit here and re-run the
% MAS build step so both stay in sync). Implements, without SimEvents and in a
% code-generation-compatible way, the message-based multi-agent black-start
% algorithm of the paper "Multi-Agent-Based Black Start of a Microgrid":
%
%   * one agent per generator (AC gens | DC gens | interface converters, in
%     buses_lf order), communication restricted to the two ring neighbours;
%   * discrete typed messages (types 0..8) with a fixed transmission delay,
%     modelled by a circular-buffer channel (one FIFO slot ring per directed
%     link, delay = round(delay/Ts) ticks);
%   * dynamic communication-graph construction (announce / follow-me /
%     I-follow-you / cannot-follow / startup-data query-reply / startup order);
%   * event-triggered secondary-control messages from each agent to its
%     follower, carrying (omega, P/Pnom, v, Q/Pnom);
%   * black-start soft start of the leader: its frequency reference ramps
%     0 -> 1 pu at framp while the voltage reference is held at 0; once the
%     frequency is nominal the voltage reference ramps 0 -> vs at vramp;
%   * agents not yet started (no connect command, breaker open) get dark
%     references (f* = 0, v* = 0): a unit awaiting its startup order must
%     not idle at nominal frequency/voltage;
%   * consensus: leader pins frequency to 1 pu; followers track the received
%     neighbour frequency (omega* channel) and share active/reactive power
%     proportionally to the nominal power MAS_Pnom (P*/Q* channels);
%   * DC agents share P like any other agent but have no frequency and no
%     reactive power: they PASS THROUGH the omega and Q data received from
%     their neighbour (both in their outgoing messages and on their inert
%     ref rows) and keep V_dc* on the v channel (vs_lf), per the framework
%     hybrid AC/DC convention;
%   * per-agent counters of secondary-control and total sent messages, plus a
%     log of the messages delivered each tick (for sequence diagrams).
%
% Inputs (elec and vmod are wired; the rest are Parameter data resolved from
% the Grid mask workspace):
%   elec         [6 x N]  rows: omega, P (pu of Sb), v (pu), Q (pu), spare,
%                         connected. For DC agents: row 1/4 echoes, row 3 V_dc.
%   vmod         [1 x Nbuses] AC bus voltage magnitudes (grid side). The
%                         PCC voltage vmod(buses_lf(i)) decides whether the
%                         communication equipment of agent i is powered
%                         (paper: comms work only above 0.5 pu at the PCC).
%   MAS_cfg_vec  [1 x 20] configuration (see layout below).
%   MAS_Pnom     [1 x N]  nominal power of each agent in pu of Sb (mBase/Sb).
%   MAS_vs       [1 x N]  per-agent voltage setpoints (= vs_lf; V_dc* for DC).
%   MAS_isdc     [1 x N]  1 for DC-generator agents, 0 otherwise.
%   MAS_buses    [1 x N]  bus index of each agent (= buses_lf).
%
% MAS_cfg_vec layout:
%    1 Ts      block sample time (s)          10 AP     trigger weight on Pn
%    2 delay   message delay (s)              11 AQ     trigger weight on Qn
%    3 leader  leader agent index             12 Tmax   max silence between SC msgs (s)
%    4 vcom    com-alive voltage threshold    13 tau    steady-state window (s), eq. (1)
%    5 kw      freq-restoration gain (1/s)    14 Assw   steady-state weight on omega
%    6 kP      P-sharing gain (1/s)           15 Assv   steady-state weight on v
%    7 kQ      Q-sharing gain (1/s)           16 Tdwell steady dwell before a query round (s)
%    8 Aw      trigger weight on omega        17 Tround min period between query rounds (s)
%    9 Av      trigger weight on v            18 vramp  black-start voltage ramp rate (pu/s)
%                                             19 Tcom   com-up dwell before announcing (s)
%                                             20 framp  leader black-start frequency ramp
%                                                       rate (pu/s; optional, default 0.2)
%
% Outputs:
%   ref        [6 x N] rows: omega*, P*, v*, Q*, connect command, connected.
%   sentSC     [1 x N] cumulative secondary-control messages sent per agent.
%   sentAll    [1 x N] cumulative messages of any type sent per agent.
%   following_o[1 x N] current parent of each agent in the graph (0 = none).
%   sending_o  [1 x N] current follower (child) of each agent (0 = none).
%   delivered  [3 x 2N] messages delivered THIS tick: [type; src; dst] per
%              column, type = -1 marks an empty column (for post-processing
%              into a message-sequence diagram).
%
% Message types (paper, Sec. "A distributed restoration procedure"):
%   0 dest>0: secondary-control data; dest==0: broadcast "new DG" announce
%   1 source wants destination to follow it
%   2 source confirms that it follows destination
%   3 source wants destination to follow it and takes the leader role
%   4 source no longer follows destination
%   5 source asks destination to stop following it / refuses a follower
%   6 source requests startup data (available capacity) from destination and
%     its neighbours (recursive query along the ring)
%   7 source sends startup data (best candidate id + capacity) to destination
%   8 startup order, routed along the ring to the optimal DG

N  = size(elec, 2);
Ts     = MAS_cfg_vec(1);  delay  = MAS_cfg_vec(2);
leader = round(MAS_cfg_vec(3));
vcom   = MAS_cfg_vec(4);
kw     = MAS_cfg_vec(5);  kP     = MAS_cfg_vec(6);  kQ = MAS_cfg_vec(7);
Aw     = MAS_cfg_vec(8);  Av     = MAS_cfg_vec(9);
AP     = MAS_cfg_vec(10); AQ     = MAS_cfg_vec(11);
Tmax   = MAS_cfg_vec(12); tau    = MAS_cfg_vec(13);
Assw   = MAS_cfg_vec(14); Assv   = MAS_cfg_vec(15);
Tdwell = MAS_cfg_vec(16); Tround = MAS_cfg_vec(17);
vramp  = MAS_cfg_vec(18); Tcom   = MAS_cfg_vec(19);
if numel(MAS_cfg_vec) >= 20
    framp = MAS_cfg_vec(20);
else
    framp = 0.2; % legacy 19-element cfg: 0 -> 1 pu in 5 s
end

D   = max(1, round(delay/Ts));  % delay in ticks
C   = 4;                        % message slots per link and tick
nL  = 2*N;                      % directed links (2 per agent: side 1 -> i-1, side 2 -> i+1)
Dss = max(1, round(tau/Ts));    % steady-state window in ticks
XI_MAX  = 0.1;                  % saturation of the frequency-restoration integrator (pu)
CHI_MAX = 2;                    % saturation of the P/Q-sharing integrators (normalised)
T_ANN   = 2;                    % re-announce period while unattached (s)
T_QUERY = 5;                    % query-round watchdog timeout (s)

persistent k chan announced tann following sending connectCmd vrampOn ...
    wrampOn wout ...
    xi chi zeta wrec Pnrec vrec Qnrec rcvValid ...
    wsent Pnsent vsent Qnsent tsent nSC nAll vout ...
    ssW ssV steadyCnt comUpCnt lastRound roundActive ...
    qActive qFrom qLeft qBestID qBestCap qBestSide qT0

if isempty(k)
    k = 0;
    chan = -ones(8, C, nL, D);      % fields: [mtype src dst p1 p2 p3 p4 p5]
    announced  = false(1, N); tann = -1e9*ones(1, N);
    comUpCnt = zeros(1, N);
    following  = zeros(1, N); sending = zeros(1, N);
    connectCmd = false(1, N); vrampOn = false(1, N);
    wrampOn = false(1, N); wout = ones(1, N);
    xi   = zeros(1, N); chi = zeros(1, N); zeta = zeros(1, N);
    wrec = ones(1, N);  Pnrec = zeros(1, N); vrec = ones(1, N); Qnrec = zeros(1, N);
    rcvValid = false(1, N);
    wsent = ones(1, N); Pnsent = zeros(1, N); vsent = ones(1, N); Qnsent = zeros(1, N);
    tsent = -1e9*ones(1, N);
    nSC = zeros(1, N); nAll = zeros(1, N);
    vout = zeros(1, N);
    ssW = zeros(Dss, N); ssV = zeros(Dss, N); steadyCnt = zeros(1, N);
    lastRound = -1e9; roundActive = false;
    qActive = false(1, N); qFrom = zeros(1, N); qLeft = zeros(1, N);
    qBestID = zeros(1, N); qBestCap = -ones(1, N); qBestSide = zeros(1, N);
    qT0 = zeros(1, N);
end

t = k*Ts;
w    = elec(1, :);
P    = elec(2, :);
v    = elec(3, :);
Q    = elec(4, :);
conn = elec(6, :) > 0.5;
Pnom = max(MAS_Pnom, 1e-9);
Pn = P./Pnom;
Qn = Q./Pnom;
isdc = MAS_isdc > 0.5;

% Communication alive: GRID-side PCC voltage above threshold (elec row 3 is
% the converter-side amplitude, ~1 pu for an initialised but disconnected
% unit, so it cannot be used; DC agents use their own V_dc since buses_lf
% only holds an AC placeholder bus for them), own generator connected, or
% being the leader (the leader has an autonomous supply).
vpcc = zeros(1, N);
for i = 1:N
    if isdc(i)
        vpcc(i) = v(i);
    else
        vpcc(i) = vmod(1, MAS_buses(i));
    end
end
com = (vpcc > vcom) | conn;
com(leader) = true;
comUpCnt = com.*(comUpCnt + 1);

% Leader startup ordered by the DSO at t = 0, energising the dead island
% with a soft start: the frequency reference ramps 0 -> 1 pu first (the
% voltage reference is held at 0 meanwhile), then the voltage ramps 0 -> vs.
if k == 0
    connectCmd(leader) = true;
    wrampOn(leader) = true;
    wout(leader) = 0;
    vrampOn(leader) = true;
    vout(leader) = 0;
end

slot = mod(k, D) + 1;

%% 1) Deliver the messages written D ticks ago into this slot, then clear it.
delivered = -ones(3, 2*N); ndel = 0;
inbox = chan(:, :, :, slot);
chan(:, :, :, slot) = -ones(8, C, nL);
for l = 1:nL
    i = floor((l-1)/2) + 1;            % sender
    s = l - 2*(i-1);                   % sender side (1 -> i-1, 2 -> i+1)
    j = ring_neigh(i, s, N);           % receiver
    r = 3 - s;                         % side of j pointing back at i
    if j == i, continue, end           % degenerate ring (N < 3)
    for c = 1:C
        mtype = inbox(1, c, l);
        if mtype < 0, continue, end
        src = round(inbox(2, c, l)); dst = round(inbox(3, c, l));
        p1 = inbox(4, c, l); p2 = inbox(5, c, l);
        p3 = inbox(6, c, l); p4 = inbox(7, c, l);
        if ~com(j), continue, end      % receiver com down: message lost
        if ndel < 2*N
            ndel = ndel + 1;
            delivered(:, ndel) = [mtype; src; dst];
        end
        inGraph = (j == leader) || (following(j) ~= 0);
        switch round(mtype)
            case 0
                if dst == 0
                    % Broadcast "new DG": offer to adopt it if I am in the
                    % graph and my follower slot is free.
                    if inGraph && sending(j) == 0 && following(j) ~= src
                        [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 1, src, 0, 0, 0, 0, 0);
                    end
                elseif dst == j && src == following(j)
                    % Secondary-control data from my parent.
                    wrec(j) = p1; Pnrec(j) = p2; vrec(j) = p3; Qnrec(j) = p4;
                    rcvValid(j) = true;
                end
            case 1 % source wants me to follow it
                if j ~= leader && following(j) == 0 && sending(j) ~= src
                    following(j) = src;
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 2, src, 0, 0, 0, 0, 0);
                end
            case 2 % source follows me
                if inGraph && following(j) ~= src && (sending(j) == 0 || sending(j) == src)
                    sending(j) = src;
                    % Bootstrap the secondary-control chain towards the new
                    % follower with the current data.
                    if isdc(j), mw = wrec(j); mq = Qnrec(j); else, mw = w(j); mq = Qn(j); end
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 0, src, mw, Pn(j), v(j), mq, 0);
                    wsent(j) = mw; Pnsent(j) = Pn(j); vsent(j) = v(j); Qnsent(j) = mq; tsent(j) = t;
                else
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 5, src, 0, 0, 0, 0, 0);
                end
            case 3 % source takes the leader role and wants me to follow it
                if following(j) ~= 0 && following(j) ~= src
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 5, src, 0, 0, 0, 0, 0);
                elseif j ~= leader
                    following(j) = src;
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 2, src, 0, 0, 0, 0, 0);
                end
            case 4 % source no longer follows me
                if sending(j) == src, sending(j) = 0; end
            case 5 % source refuses me / asks me to stop following it
                if following(j) == src
                    following(j) = 0; rcvValid(j) = false;
                    if sending(j) ~= 0
                        % Warn my own follower: the subtree dissolves.
                        sside = side_of(j, sending(j), N);
                        [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, sside, 5, sending(j), 0, 0, 0, 0, 0);
                        sending(j) = 0;
                    end
                end
            case 6 % startup-data request
                [cid, cap] = local_candidate(j, com(j), conn(j), connectCmd(j), following(j), Pnom(j));
                if qActive(j)
                    % Busy with another query: answer with local data only.
                    [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 7, src, cid, cap, 0, 0, 0);
                else
                    qActive(j) = true; qFrom(j) = r; qT0(j) = t;
                    qBestID(j) = cid; qBestCap(j) = cap; qBestSide(j) = 0;
                    cnt = 0;
                    for side = 1:2
                        nb = ring_neigh(j, side, N);
                        if side ~= r && nb ~= j && nb ~= src && com(nb)
                            [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, side, 6, nb, 0, 0, 0, 0, 0);
                            cnt = cnt + 1;
                        end
                    end
                    qLeft(j) = cnt;
                    if cnt == 0
                        [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, r, 7, src, qBestID(j), qBestCap(j), 0, 0, 0);
                        qActive(j) = false;
                    end
                end
            case 7 % startup-data reply: p1 = candidate id, p2 = capacity
                if qActive(j) && qLeft(j) > 0
                    qLeft(j) = qLeft(j) - 1;
                    if p1 > 0 && p2 > qBestCap(j)
                        qBestID(j) = round(p1); qBestCap(j) = p2; qBestSide(j) = r;
                    end
                    if qLeft(j) == 0
                        if qFrom(j) > 0 % relay the aggregate answer upstream
                            [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, qFrom(j), 7, ring_neigh(j, qFrom(j), N), qBestID(j), qBestCap(j), 0, 0, 0);
                            qActive(j) = false;
                        else            % query originator (the leader)
                            qActive(j) = false; roundActive = false;
                            if qBestID(j) > 0 && qBestID(j) ~= j && qBestSide(j) > 0
                                [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, qBestSide(j), 8, ring_neigh(j, qBestSide(j), N), qBestID(j), N, 0, 0, 0);
                            end
                        end
                    end
                end
            case 8 % startup order: p1 = target agent, p2 = remaining hops
                if round(p1) == j
                    connectCmd(j) = true;
                    vrampOn(j) = vpcc(j) < 0.3; % energising a dead grid: ramp v*
                    if vrampOn(j), vout(j) = max(vpcc(j), 0); end
                elseif p2 > 0
                    fside = 3 - r; % keep travelling around the ring
                    nb = ring_neigh(j, fside, N);
                    if nb ~= j && com(nb)
                        [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, j, fside, 8, nb, p1, p2 - 1, 0, 0, 0);
                    end
                end
            otherwise
                % unknown type: ignore
        end
    end
end

%% 2) Announce new agents (broadcast type 0) while unattached, with retry.
% The com equipment must have been powered for Tcom seconds before the first
% announce (this also rejects the phantom initial bus voltages of the
% load-flow-initialised network at t = 0).
for i = 1:N
    if i ~= leader && com(i) && comUpCnt(i)*Ts >= Tcom ...
            && following(i) == 0 && (t - tann(i)) >= T_ANN
        sent_any = false;
        for side = 1:2
            nb = ring_neigh(i, side, N);
            if nb ~= i && com(nb)
                [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, i, side, 0, 0, 0, 0, 0, 0, 0);
                sent_any = true;
            end
        end
        if sent_any
            announced(i) = true; tann(i) = t;
        end
    end
end

%% 3) Steady-state detection, eq. (1): psi = Assv*(v(t)-v(t-tau))^2 + ...
idxss = mod(k, Dss) + 1;
wOld = ssW(idxss, :); vOld = ssV(idxss, :);
ssW(idxss, :) = w; ssV(idxss, :) = v;
if k >= Dss
    psi = Assw*(w - wOld).^2 + Assv*(v - vOld).^2;
else
    psi = 1e9*ones(1, N);
end
steadyCnt = (psi < 1).*(steadyCnt + 1);

%% 4) Leader query rounds (types 6/7/8): one after every reached steady state,
%% plus a periodic retry every Tround while steady.
dwellTicks = max(1, round(Tdwell/Ts));
% Watchdog: a reply was lost (com went down mid-query). Close the pending
% query with the best data known so far, so no agent stays blocked.
for i = 1:N
    if qActive(i) && (t - qT0(i)) > T_QUERY
        if qFrom(i) > 0
            [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, i, qFrom(i), 7, ring_neigh(i, qFrom(i), N), qBestID(i), qBestCap(i), 0, 0, 0);
        elseif i == leader
            roundActive = false;
        end
        qActive(i) = false;
    end
end
if conn(leader) && ~roundActive && ~qActive(leader) && steadyCnt(leader) >= dwellTicks ...
        && (steadyCnt(leader) == dwellTicks || (t - lastRound) >= Tround)
    qActive(leader) = true; qFrom(leader) = 0; qT0(leader) = t;
    qBestID(leader) = 0; qBestCap(leader) = -1; qBestSide(leader) = 0;
    cnt = 0;
    for side = 1:2
        nb = ring_neigh(leader, side, N);
        if nb ~= leader && com(nb)
            [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, leader, side, 6, nb, 0, 0, 0, 0, 0);
            cnt = cnt + 1;
        end
    end
    qLeft(leader) = cnt; lastRound = t;
    if cnt == 0
        qActive(leader) = false;
    else
        roundActive = true;
    end
end

%% 5) Event-triggered secondary-control messages to the follower.
for i = 1:N
    if conn(i) && sending(i) ~= 0 && com(i)
        if isdc(i)
            mw = wrec(i); mq = Qnrec(i); % DC: pass through neighbour f and Q
        else
            mw = w(i); mq = Qn(i);
        end
        trig = Aw*(mw - wsent(i))^2 + Av*(v(i) - vsent(i))^2 ...
             + AP*(Pn(i) - Pnsent(i))^2 + AQ*(mq - Qnsent(i))^2;
        if trig > 1 || (t - tsent(i)) > Tmax
            side = side_of(i, sending(i), N);
            [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, i, side, 0, sending(i), mw, Pn(i), v(i), mq, 0);
            wsent(i) = mw; Pnsent(i) = Pn(i); vsent(i) = v(i); Qnsent(i) = mq;
            tsent(i) = t;
        end
    end
end

%% 6) Consensus integrators (Euler, dt = Ts).
for i = 1:N
    if conn(i)
        if i == leader
            if ~wrampOn(i) % frozen during the black-start frequency ramp
                xi(i) = sat(xi(i) + Ts*kw*(1 - w(i)), XI_MAX);
            end
        elseif rcvValid(i)
            chi(i) = sat(chi(i) + Ts*kP*(Pnrec(i) - Pn(i)), CHI_MAX);
            if ~isdc(i)
                xi(i)   = sat(xi(i)   + Ts*kw*(wrec(i) - w(i)),  XI_MAX);
                zeta(i) = sat(zeta(i) + Ts*kQ*(Qnrec(i) - Qn(i)), CHI_MAX);
            end
        end
    end
end

%% 7) Frequency and voltage references (black-start soft start: f, then v).
for i = 1:N
    if wrampOn(i)
        wout(i) = min(1, wout(i) + Ts*framp);
        if wout(i) >= 1, wrampOn(i) = false; end
    end
    if ~connectCmd(i) && ~conn(i)
        vout(i) = 0; % dark until the startup order arrives
    elseif vrampOn(i)
        if ~wrampOn(i) % voltage held at 0 until the frequency is nominal
            vout(i) = min(MAS_vs(i), vout(i) + Ts*vramp);
            if vout(i) >= MAS_vs(i), vrampOn(i) = false; end
        end
    else
        vout(i) = MAS_vs(i);
    end
end

%% 8) Outputs.
ref = zeros(6, N);
for i = 1:N
    if isdc(i)
        if rcvValid(i), ref(1, i) = wrec(i); else, ref(1, i) = 1; end
        ref(4, i) = 0;                 % Q channel inert for DC agents
    else
        if ~connectCmd(i) && ~conn(i)
            ref(1, i) = 0; % dark until the startup order arrives
        elseif wrampOn(i)
            ref(1, i) = wout(i);
        else
            ref(1, i) = 1 + xi(i);
        end
        ref(4, i) = zeta(i)*Pnom(i);
    end
    ref(2, i) = chi(i)*Pnom(i);
    ref(3, i) = vout(i);
    ref(5, i) = double(connectCmd(i));
    ref(6, i) = double(conn(i));
end
sentSC = nSC;
sentAll = nAll;
following_o = following;
sending_o = sending;
k = k + 1;
end

function y = sat(x, lim)
% Symmetric saturation.
y = min(max(x, -lim), lim);
end

function nb = ring_neigh(i, side, N)
% Ring neighbour of agent i: side 1 -> i-1, side 2 -> i+1 (with wrap-around).
if side == 1
    nb = i - 1;
    if nb < 1, nb = N; end
else
    nb = i + 1;
    if nb > N, nb = 1; end
end
end

function side = side_of(i, target, N)
% Side of agent i on which ring neighbour `target` sits (0 if not adjacent).
if target == ring_neigh(i, 1, N)
    side = 1;
elseif target == ring_neigh(i, 2, N)
    side = 2;
else
    side = 0;
end
end

function [cid, cap] = local_candidate(j, comj, connj, cmdj, followingj, Pnomj)
% Startup candidate reported by agent j: itself, if it is communicating, in
% the graph, and not yet (commanded to be) connected. Capacity = nominal power.
if comj && ~connj && ~cmdj && followingj ~= 0
    cid = j; cap = Pnomj;
else
    cid = 0; cap = -1;
end
end

function [chan, nAll, nSC] = do_send(chan, nAll, nSC, slot, C, N, com, i, side, mtype, dst, p1, p2, p3, p4, p5)
% Queue a message from agent i on the given ring side. It will be delivered
% after the configured delay (the write slot is read D ticks later). Messages
% are dropped if the sender is com-down, the side is invalid, or the per-tick
% link buffer is full.
if side < 1 || side > 2 || ~com(i)
    return
end
l = 2*(i - 1) + side;
for c = 1:C
    if chan(1, c, l, slot) < 0
        chan(:, c, l, slot) = [mtype; i; dst; p1; p2; p3; p4; p5];
        nAll(i) = nAll(i) + 1;
        if mtype == 0 && dst > 0
            nSC(i) = nSC(i) + 1;
        end
        return
    end
end
end
