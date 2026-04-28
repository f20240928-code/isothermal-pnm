%% =========================================================================
%%  dryingLoop_3D.m
%%  Exact 3D physics replica of the prof's main drying loop.
%%
%%  PHYSICS – NOTHING CHANGED.  Every formula is identical to the 2D code.
%%  The only adaptations are:
%%    (a) Index ranges extended to 3D counts  (Nt, NT, NP, Nsurf, m, n, p)
%%    (b) Vtslice sum covers Y + X + Z throats  (3D slice volume)
%%    (c) Call labelLiquidClusters_3D  instead of dry_e_hoshen_kopelmann
%%    (d) network.dynamic fields updated each step for video
%%    (e) Video writer uses a 3D visualisation call
%%
%%  REQUIRES (workspace – set by main script before calling this file)
%%  ------------------------------------------------------------------
%%  From network.static (unpacked at top):
%%    NP, NT, Nt, Nty, Ntx, Ntz, m, n, p, mbl, Nsurf
%%    pnp, pnt, tnp, tnt, tnp2
%%    rt, At, Lt, Vt, Vtslice, TotalVolume_Pore_Network
%%  Physical constants (set in main / params file):
%%    Patm, sigma, mu, RhoL, Peqv
%%  Conductances (set before loop):
%%    gv(NT,1), ginf(NP,1)
%%  State arrays (initialised before loop):
%%    ts(NT,1)   – throat saturation level  (1=full liquid, 0=vapour)
%%    ps(NP,1)   – pore state               (1=liquid, 0=vapour)
%%  From solveVaporPressure_3D  (must be run once before this loop):
%%    A0, B0, Totaltime, Tot_CT
%%    Atactive, rcandidate, Pressure_meniscus, gw0, Pore1, Pore2
%%  From solveLiquidPressure_3D  (must be run once before this loop):
%%    Peq, A0_liquid, B0_liquid, Volume_Liquid, Total_Cross_A
%%    NumberCluster, Number_Cluster, q, Sslice
%%  From identifyCandidateThroats_3D  (must be run once before this loop):
%%    MM, Pore_Cond
%%  From labelLiquidClusters_3D  (must be run once before this loop):
%%    tl(Nt,1)
%%  Video:
%%    Obj  – VideoWriter object, already open
%%
%%  PRODUCES
%%  --------
%%  Mev(Nt-1,1)              : total evaporation rate per step  [kg/s]
%%  TotalVolume_Liquid(Nt-1) : total liquid volume per step     [m³]
%%  Ss(Nt-1,1)               : global saturation per step       [0-1]
%%  timee(Nt-1,1)            : cumulative simulation time       [s]
%%  Sslice(101, m-1)         : saturation profiles per Y-slice
%%  Comptime(Nt-1,1)         : wall-clock time per step         [s]
%%  network.dynamic          : updated every step for video
%% =========================================================================

%% --- Unpack static fields once (avoids repeated struct access in loop) ---
NP     = network.static.NP;
NT     = network.static.NT;
Nt     = network.static.Nt;
Nty    = network.static.Nty;
Ntx    = network.static.Ntx;
Ntz    = network.static.Ntz;
m      = network.static.m;
n      = network.static.n;
p      = network.static.p;
Nsurf  = network.static.Nsurf;
pnp    = network.static.pnp;    % NP × 6
pnt    = network.static.pnt;    % NP × 6
tnp    = network.static.tnp;    % NT × 2
At     = network.static.At;     % NT × 1
rt     = network.static.rt;     % Nt × 1
Lt     = network.static.Lt;     % NT × 1
Vt     = network.static.Vt;     % Nt × 1
Vtslice                = network.static.Vtslice;          % 1 × (m-1)
TotalVolume_Pore_Network = network.static.TotalVolume;

%% --- Pre-allocate output tracking arrays --------------------------------
Mev                = zeros(NT, 1);    % momentary evap rate (indexed by step s)
TotalVolume_Liquid = zeros(Nt-1, 1);
Ss                 = zeros(Nt-1, 1);
timee              = zeros(Nt-1, 1);
Comptime           = zeros(Nt-1, 1);
Ev                 = zeros(NP, 1);    % per-pore evaporation rate

%% =========================================================================
%%  MAIN DRYING LOOP
%%  2D: for s = 1:Nt-1   →   3D: identical range (Nt is now 3D interior count)
%% =========================================================================
for s = 1:Nt-1

    P_liquid = zeros(NP, 1);

    %======================================================================
    %  STEP 1 – Solve vapour pressures  (sparse linear system Ax = b)
    %  Identical to 2D.  index picks only unknown vapour pores (ps==0, Peq==0).
    %  index_liquid picks all liquid pores.
    %======================================================================
    index         = find(ps == 0 & Peq == 0);
    index_liquid  = find(ps == 1);
    NPw           = length(index_liquid);

    % --- Liquid pressure sub-system ---------------------------------------
    % 2D: if NPw ~= 0 ... A_liquid\B_liquid ... end   UNCHANGED.
    if NPw ~= 0
        A_liquid = A0_liquid(index_liquid, index_liquid);
        B_liquid = B0_liquid(index_liquid);
        X_liquid = A_liquid \ B_liquid;
        P_liquid(index_liquid) = X_liquid;
    end

    % --- Vapour pressure sub-system ---------------------------------------
    % 2D: A = A0(index,index); B = B0(index); X = A\B   UNCHANGED.
    A = A0(index, index);
    B = B0(index);
    X = A \ B;

    % Convert log-activity solution back to partial pressure
    % 2D: P_X_Vapour = (1 - exp(X)) * Patm   UNCHANGED.
    P_X_Vapour = (1 - exp(X)) * Patm;

    % Assemble full vapour pressure vector
    % 2D: Pv = Peq;  Pv(index) = P_X_Vapour   UNCHANGED.
    Pv        = Peq;
    Pv(index) = P_X_Vapour;

    %======================================================================
    %  STEP 2 – Evaporation rate at the top boundary layer surface
    %  2D: for i = NP-Nsurf+1 : NP   [Nsurf = n in 2D]
    %  3D: Nsurf = n*p  (top face is a 2D slab)   range is IDENTICAL in form.
    %======================================================================
    for i = NP - Nsurf + 1 : NP
        % 2D: Ev(i) = -ginf(i)*Patm*log((Patm - Pv(i))/Patm)   UNCHANGED.
        Ev(i) = -ginf(i) * Patm * log((Patm - Pv(i)) / Patm);
    end

    % 2D: Mev(s) = sum(Ev)   UNCHANGED.
    Mev(s) = sum(Ev);

    %======================================================================
    %  STEP 3 – Vapour diffusive flow at every throat
    %  2D: Flow(i) = gv(i)*Patm*abs(log(...))   UNCHANGED.
    %======================================================================
    Flow = zeros(NT, 1);
    for i = 1:NT
        Flow(i) = gv(i) * Patm * ...
            abs(log((Patm - Pv(tnp(i,1))) / (Patm - Pv(tnp(i,2)))));
    end

    %======================================================================
    %  STEP 4 – Distribute flow to throats via active cross-section weighting
    %  Identical to 2D in formula.  Pore1, Pore2, Total_Cross_A, Atactive
    %  are all NT-length arrays that already account for 3D topology.
    %======================================================================
    Total_Flows_Pore1 = zeros(NT, 1);
    Total_Flows_Pore2 = zeros(NT, 1);
    Total_Mev_Throat  = zeros(NT, 1);

    for i = 1:NT
        % Contribution via Pore1
        % 2D: if Total_Cross_A(Pore1(i)) ~= 0 ...   UNCHANGED.
        if Total_Cross_A(Pore1(i)) ~= 0
            Total_Flows_Pore1(i) = ...
                (Atactive(i) / Total_Cross_A(Pore1(i))) * ...
                sum(Flow(nonzeros(pnt(Pore1(i), :))));
        end

        % Contribution via Pore2
        if Total_Cross_A(Pore2(i)) ~= 0
            Total_Flows_Pore2(i) = ...
                (Atactive(i) / Total_Cross_A(Pore2(i))) * ...
                sum(Flow(nonzeros(pnt(Pore2(i), :))));
        end

        % 2D: Total_Mev_Throat(i) = Total_Flows_Pore1(i) + Total_Flows_Pore2(i)
        Total_Mev_Throat(i) = Total_Flows_Pore1(i) + Total_Flows_Pore2(i);
    end

    % 2D: TMT = Total_Mev_Throat(find(Total_Mev_Throat ~= 0))   UNCHANGED.
    TMT = Total_Mev_Throat(Total_Mev_Throat ~= 0);

    %======================================================================
    %  STEP 5 – Per-cluster drying time and maximum-radius throat selection
    %  Identical to 2D.  tl, NumberCluster, rcandidate, gw0 are all
    %  Nt-length arrays consistent with 3D interior-throat indexing.
    %======================================================================
    Mevv             = zeros(NumberCluster, 1);
    max_radius_cluster = zeros(NumberCluster, 1);
    DTT              = zeros(NumberCluster, 1);
    indexmaxthroats  = zeros(NumberCluster, 1);
    Mass_Flow_liquid = zeros(NumberCluster, 1);

    for pp = 1:NumberCluster

        % Total evaporation rate from cluster pp
        % 2D: Mevv(p) = sum(Total_Mev_Throat(tl==p))   UNCHANGED.
        Mevv(pp) = sum(Total_Mev_Throat(tl == pp));

        % Throats belonging to cluster pp
        indexthroats = find(tl == pp);

        % Throat with largest candidate radius in cluster
        % 2D: [max_radius_cluster(p), pos1] = max(rcandidate(indexthroats))   UNCHANGED.
        [max_radius_cluster(pp), pos1] = max(rcandidate(indexthroats));
        indexmaxthroats(pp) = indexthroats(pos1);

        % Liquid mass flow to the meniscus throat
        % 2D: two-branch if/elseif on ps of Pore1/Pore2   UNCHANGED.
        imt = indexmaxthroats(pp);
        if (ps(Pore1(imt)) == 0 && ps(Pore2(imt)) == 1)
            Mass_Flow_liquid(pp) = gw0(imt) * ...
                abs(Pressure_meniscus(imt) - P_liquid(tnp(imt, 2)));
        elseif (ps(Pore1(imt)) == 1 && ps(Pore2(imt)) == 0)
            Mass_Flow_liquid(pp) = gw0(imt) * ...
                abs(Pressure_meniscus(imt) - P_liquid(tnp(imt, 1)));
        end

        % Physical consistency check (same as 2D)
        if Mass_Flow_liquid(pp) + Total_Mev_Throat(imt) < Mevv(pp)
            error('[dryingLoop_3D] Step %d cluster %d: liquid flow < evaporation!', s, pp);
        end

        % Drying time for this cluster
        % 2D: if Mevv(p) ~= 0; DTT(p) = (Volume_Liquid(imt)*RhoL)/Mevv(p); end
        if Mevv(pp) ~= 0
            DTT(pp) = (Volume_Liquid(indexmaxthroats(pp)) * RhoL) / Mevv(pp);
        end
    end

    %======================================================================
    %  STEP 6 – Singleton (non-cluster) candidate throats and their times
    %  2D: NN = find(tl(MM)==0 & Total_Mev_Throat(MM)~=0)   UNCHANGED.
    %======================================================================
    NN = find(tl(MM) == 0 & Total_Mev_Throat(MM) ~= 0);

    Tt = zeros(length(NN), 1);
    for i = 1:length(NN)
        % 2D: Tt(i) = (Volume_Liquid(MM(NN(i)))*RhoL)/Total_Mev_Throat(MM(NN(i)))
        Tt(i) = (Volume_Liquid(MM(NN(i))) * RhoL) / ...
                 Total_Mev_Throat(MM(NN(i)));
    end

    %======================================================================
    %  STEP 7 – Global minimum drying time dt
    %  2D: Candidate_min, T_tot, [dt indices] = min(nonzeros(T_tot))   UNCHANGED.
    %======================================================================
    TIT            = find(DTT ~= 0);
    Candidate_min  = [MM(NN); indexmaxthroats(TIT)];
    T_tot          = [Tt; DTT(TIT)];

    [dt, ~] = min(nonzeros(T_tot));

    %======================================================================
    %  STEP 8 – Global saturation tracking
    %  2D: TotalVolume_Liquid(s) = sum(Volume_Liquid)   UNCHANGED.
    %======================================================================
    TotalVolume_Liquid(s) = sum(Volume_Liquid);
    Ss(s) = TotalVolume_Liquid(s) / TotalVolume_Pore_Network;

    %======================================================================
    %  STEP 9 – Saturation slice profiles  (saved every 1 % saturation)
    %
    %  2D formula (prof):
    %    Vliquidslice =
    %      sum(reshape(Vt(1:n*(m-1)) .* ts(1:n*(m-1)),  n, m-1), 1)
    %    + sum(reshape(Vt(n*(m-1)+1:n*(m-1)+(n-1)*(m-1)) .*
    %                  ts(n*(m-1)+1:n*(m-1)+(n-1)*(m-1)), n-1, m-1), 1)
    %
    %  3D extension: add Z-throat contribution; all three directions
    %  use the same per-slice index layout built in the network generator.
    %
    %  Nty_per_slice = n*p,  Ntx_per_slice = (n-1)*p,  Ntz_per_slice = n*(p-1)
    %======================================================================
    if floor(100 * Ss(s)) < q

        Nty_per_slice = n * p;
        Ntx_per_slice = (n-1) * p;
        Ntz_per_slice = n * (p-1);

        Vliquidslice = zeros(1, m-1);
        for sl = 1:m-1
            % Y-throats for this slice
            yIdx = (sl-1)*Nty_per_slice + 1 : sl*Nty_per_slice;
            % X-throats for this slice
            xIdx = Nty + (sl-1)*Ntx_per_slice + 1 : Nty + sl*Ntx_per_slice;
            % Z-throats for this slice  (new in 3D)
            zIdx = Nty+Ntx + (sl-1)*Ntz_per_slice + 1 : Nty+Ntx + sl*Ntz_per_slice;

            Vliquidslice(sl) = sum(Vt(yIdx) .* ts(yIdx)) ...
                             + sum(Vt(xIdx) .* ts(xIdx)) ...
                             + sum(Vt(zIdx) .* ts(zIdx));
        end

        % 2D: Sslice(100-q+1,:) = Vliquidslice ./ Vtslice   UNCHANGED.
        Sslice(100 - q + 1, :) = Vliquidslice ./ Vtslice;
        q = q - 1;
        fprintf('  Saturation profile saved at Ss = %.3f  (q=%d)\n', Ss(s), q);
    end

    fprintf('[dryingLoop_3D] s=%d  dt=%.4e s  Ss=%.4f  Clusters=%d\n', ...
            s, dt, Ss(s), NumberCluster);

    %======================================================================
    %  STEP 10 – Advance time
    %  2D: Totaltime = Totaltime + dt;  timee(s) = Totaltime   UNCHANGED.
    %======================================================================
    Totaltime  = Totaltime + dt;
    timee(s)   = Totaltime;

    %======================================================================
    %  STEP 11 – Update throat states for singleton candidates
    %  2D formula for ts update:
    %    ts(t) = ts(t) - (dt * TMT(t) / RhoL) / Vt(t)
    %  then zero-out or update dependent arrays.   UNCHANGED.
    %======================================================================
    for i = 1:length(NN)
        t = MM(NN(i));   % throat index

        ts(t) = ts(t) - ((dt * Total_Mev_Throat(t) / RhoL) / Vt(t));

        if ts(t) < 1e-15
            ts(t) = 0;
        end

        if ts(t) == 0
            % Throat fully evaporated
            Atactive(t)          = 0;
            Pressure_meniscus(t) = 0;
            Volume_Liquid(t)     = 0;
            gw0(t)               = 0;
            rcandidate(t)        = 0;
        else
            % Throat partially filled
            Atactive(t)      = At(t);
            Volume_Liquid(t) = ts(t) * Vt(t);
        end
    end

    %======================================================================
    %  STEP 12 – Update throat states for cluster max-radius throats
    %  Same formula as STEP 11 but uses Mevv(p) (cluster evap rate).
    %  Adds gw0 update for partially filled throats.   UNCHANGED.
    %======================================================================
    for pp = 1:NumberCluster
        t = indexmaxthroats(pp);

        ts(t) = ts(t) - ((dt * Mevv(pp) / RhoL) / Vt(t));

        if ts(t) < 1e-15
            ts(t) = 0;
        end

        if ts(t) == 0
            % Throat fully evaporated
            Atactive(t)          = 0;
            Pressure_meniscus(t) = 0;
            Volume_Liquid(t)     = 0;
            gw0(t)               = 0;
            rcandidate(t)        = 0;
        else
            % Throat partially filled  (gw0 also updated here, same as 2D)
            Atactive(t)      = At(t);
            Volume_Liquid(t) = ts(t) * Vt(t);
            gw0(t)           = ((pi * rt(t).^4) ./ ...
                                (mu * 8 * Lt(t) * ts(t))) * RhoL;
        end
    end

    %======================================================================
    %  STEP 13 – Cluster relabelling (replaces dry_e_hoshen_kopelmann)
    %  labelLiquidClusters_3D is state-driven: call it unconditionally
    %  after every ts update.  It refreshes tl, NumberCluster,
    %  ps, Tot_CT, Peq, A0_liquid, B0_liquid, B0, Total_Cross_A.
    %======================================================================
    computertime2 = clock;

    labelLiquidClusters_3D;

    Number_Cluster(s) = NumberCluster;
    Comptime(s)       = etime(clock, computertime2);

    %======================================================================
    %  STEP 14 – Refresh candidate list
    %  2D: MM = find(rcandidate ~= 0)   UNCHANGED.
    %======================================================================
    MM = find(rcandidate ~= 0);

    %======================================================================
    %  STEP 15 – Update network.dynamic and write video frame
    %  (replaces prof's Networkdrawing4nonperiodic + writeVideo)
    %======================================================================
    network.dynamic.Pressure(ps == 1)  = P_liquid(ps == 1);
    network.dynamic.Pressure(ps == 0)  = Pv(ps == 0);
    network.dynamic.Saturation         = ps;
    network.dynamic.FlowRate(1:NT)     = Flow;
    network.dynamic.Velocity(1:Nt)     = Flow(1:Nt) ./ max(At(1:Nt), eps);
    network.dynamic.Occupancy(1:Nt)    = double(ts(1:Nt) > 0);
    network.dynamic.Time               = Totaltime;
    network.dynamic.Step               = s;

    % Write video frame using 3D network drawing
    networkDrawing3D;
    F = getframe(gcf);
    writeVideo(Obj, F);
    clf;

end  % end main drying loop

close(Obj);
fprintf('[dryingLoop_3D]  Simulation complete.  Steps: %d  FinalSs: %.4f\n', ...
        Nt-1, Ss(Nt-1));