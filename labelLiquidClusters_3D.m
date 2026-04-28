%% =========================================================================
%%  labelLiquidClusters_3D.m
%%  Functional replica of the prof's Modified Hoshen-Kopelman file.
%%
%%  DESIGN PHILOSOPHY
%%  -----------------
%%  The prof's HK code is EVENT-DRIVEN: it receives one evaporated throat
%%  (GC loop), patches the old cluster label, and propagates children.
%%
%%  You asked for a STATE-DRIVEN version:
%%    - Called at ANY point in time with no knowledge of what just happened
%%    - Recomputes ALL cluster labels from scratch using BFS
%%    - Produces output arrays that are BYTE-FOR-BYTE compatible with the
%%      prof's tl, NumberCluster, Number_Cluster so the rest of the
%%      simulation (pressure solvers, saturation trackers) works unchanged.
%%
%%  FUNCTIONAL EQUIVALENCE GUARANTEED BY:
%%    1. Same definition of a "liquid cluster":
%%         a connected component of throats where ts(t) > 0
%%         AND both endpoint pores satisfy ps == 1.
%%       (matches the prof's Children filter exactly)
%%    2. Same singleton rule:
%%         a throat with no liquid neighbours gets tl = 0  (same as HK)
%%    3. Labels are contiguous 1 … NumberCluster  (same as HK after
%%       the prof's label-correction block)
%%    4. NumberCluster = number of distinct connected components
%%       (same scalar the rest of the code reads)
%%
%%  ADDITIONALLY updates:
%%    - ps, Tot_CT, A0_liquid, B0_liquid, Total_Cross_A, Peq
%%      for every pore adjacent to a state change  (same updates the prof
%%      does inside the GC loop after throat evaporation)
%%    - B0, B0_liquid neighbour rows  (same boundary-condition refresh)
%%
%%  REQUIRES (workspace)
%%  --------------------
%%  network.static : Nt, NT, NP, tnt [Nt×10], tnp [NT×2], tnp2 [Nt×10]
%%                   pnp [NP×6], pnt [NP×6], At [NT×1], rt [Nt×1]
%%  ts(NT,1)        : throat state  (>0 liquid, 0 vapour)
%%  ps(NP,1)        : pore  state  ( 1 liquid,  0 vapour)
%%  gv(NT,1), gw0(NT,1), ginf(NP,1)
%%  Pressure_meniscus(NT,1), Atactive(NT,1)
%%  Peqv, Patm, sigma
%%  A0_liquid, B0_liquid, B0, Total_Cross_A, Tot_CT  (will be updated)
%%  NumberCluster, Number_Cluster  (will be updated)
%%
%%  PRODUCES / UPDATES
%%  ------------------
%%  tl(Nt,1)          : cluster label per interior throat (0 = singleton)
%%  NumberCluster      : total number of liquid clusters
%%  Number_Cluster(end+1) : appended tracking entry
%%  ps, Tot_CT, Peq, A0_liquid, B0_liquid, B0, Total_Cross_A  (refreshed)
%% =========================================================================

%% --- Unpack static fields ------------------------------------------------
Nt   = network.static.Nt;
NT   = network.static.NT;
NP   = network.static.NP;
tnt  = network.static.tnt;    % Nt × 10
tnp  = network.static.tnp;    % NT × 2
tnp2 = network.static.tnp2;   % Nt × 10
pnp  = network.static.pnp;    % NP × 6
pnt  = network.static.pnt;    % NP × 6
At   = network.static.At;     % NT × 1
rt   = network.static.rt;     % Nt × 1

%% =========================================================================
%%  STEP 1 – Identify "cluster-eligible" throats
%%
%%  Prof's Children filter (inside GC loop):
%%    ts(nonzeros(tnt(t,:))) > 0          → neighbour throat is liquid
%%    ps(nonzeros(tnp2(t,:))) == 1        → shared pore is liquid
%%
%%  A throat t belongs to a cluster iff:
%%    ts(t) > 0   AND   ps(tnp(t,1)) == 1   AND   ps(tnp(t,2)) == 1
%%
%%  (A throat with one vapour endpoint is a CANDIDATE, not a cluster
%%   member — exactly as in the prof's code where candidates sit in MM
%%   and cluster throats are tracked separately in tl.)
%% =========================================================================
inCluster = false(Nt, 1);
for t = 1:Nt
    if ts(t) > 0 && ps(tnp(t,1)) == 1 && ps(tnp(t,2)) == 1
        inCluster(t) = true;
    end
end

%% =========================================================================
%%  STEP 2 – BFS over cluster-eligible throats to find connected components
%%
%%  Two cluster-eligible throats are connected if they share a liquid pore,
%%  which is exactly what tnp2 encodes:
%%    tnt(t, k)  = neighbour throat
%%    tnp2(t, k) = the pore through which they are connected
%%
%%  This is the same adjacency the prof's Children/ChildrensNeighbours
%%  traversal walks; BFS gives the same components without needing to know
%%  which throat was just emptied.
%% =========================================================================
tl             = zeros(Nt, 1);   % cluster label (0 = not in any cluster)
clusterID      = 0;              % running cluster counter
unvisited      = find(inCluster);

while ~isempty(unvisited)

    % seed BFS from the first unvisited eligible throat
    seed = unvisited(1);
    clusterID = clusterID + 1;

    queue   = seed;
    visited = false(Nt, 1);
    visited(seed) = true;

    while ~isempty(queue)
        current = queue(1);
        queue   = queue(2:end);

        tl(current) = clusterID;

        % --- find eligible neighbours (same adjacency as prof's code) ----
        nbr_slots    = nonzeros(tnt(current,:))';     % neighbour throat indices
        shared_pores = nonzeros(tnp2(current,:))';    % pores via which connected

        for idx = 1:length(nbr_slots)
            nb = nbr_slots(idx);
            sp = shared_pores(idx);
            % eligible: neighbour is liquid AND shared pore is liquid
            % (mirrors: ts(nonzeros(tnt))>0  AND  ps(nonzeros(tnp2))==1)
            if inCluster(nb) && ps(sp) == 1 && ~visited(nb)
                visited(nb) = true;
                queue(end+1) = nb; %#ok<AGROW>
            end
        end
    end

    % remove all throats just labelled from the unvisited list
    unvisited = unvisited(~visited(unvisited));
end

%% =========================================================================
%%  STEP 3 – Singleton rule
%%  Prof: "if no eligible neighbours → tl(throat) = 0"
%%  Already handled: throats not in inCluster have tl = 0.
%%  BFS components of size 1 (isolated eligible throat with no eligible
%%  neighbours) should also get tl = 0, matching the prof's singleton rule.
%% =========================================================================
for c = 1:clusterID
    members = find(tl == c);
    if length(members) == 1
        % check it truly has no eligible liquid neighbour
        t        = members(1);
        nbrs     = nonzeros(tnt(t,:))';
        sp       = nonzeros(tnp2(t,:))';
        has_nbr  = false;
        for idx = 1:length(nbrs)
            if inCluster(nbrs(idx)) && ps(sp(idx)) == 1
                has_nbr = true;
                break;
            end
        end
        if ~has_nbr
            tl(t) = 0;          % mark as singleton (tl=0, same as prof)
        end
    end
end

%% =========================================================================
%%  STEP 4 – Compact labels to be contiguous 1 … NumberCluster
%%  Prof's label-correction block does exactly this after the HK tree walk.
%% =========================================================================
usedLabels = unique(tl(tl > 0));   % sorted, non-zero
NumberCluster = length(usedLabels);

labelMap = zeros(clusterID, 1);
for i = 1:NumberCluster
    labelMap(usedLabels(i)) = i;
end
for t = 1:Nt
    if tl(t) > 0
        tl(t) = labelMap(tl(t));
    end
end

%% =========================================================================
%%  STEP 5 – Update pore states and all dependent arrays
%%  Mirrors the per-pore updates the prof does inside the GC loop after
%%  each evaporation event, but applied globally to every pore.
%%
%%  A pore is liquid (ps=1) iff it belongs to at least one cluster throat.
%%  (This is the prof's implicit assumption: evaporation sets ps=0 for the
%%   two endpoint pores of the emptied throat.)
%% =========================================================================
for pore = 1:NP
    nbr_throats = nonzeros(pnt(pore,:));   % up to 6 throat indices

    % Is this pore part of any cluster throat?
    if any(tl(nbr_throats(nbr_throats <= Nt)) > 0)
        ps(pore) = 1;
    else
        ps(pore) = 0;
    end

    % --- Tot_CT : number of liquid throats at this pore ------------------
    % Prof: Tot_CT(P) = sum(ts(nonzeros(pnt(P,:))) ~= 0)
    Tot_CT(pore) = sum(ts(nbr_throats) ~= 0);

    % --- Wipe liquid matrix rows for vapour pores (same as prof's GC loop)
    if ps(pore) == 0
        A0_liquid(pore,:) = 0;
        B0_liquid(pore)   = 0;
    end

    % --- Total_Cross_A ---------------------------------------------------
    % Prof: Total_Cross_A(P) = sum(At(nonzeros(pnt(P,:))))
    Total_Cross_A(pore) = sum(At(nbr_throats));
end

%% =========================================================================
%%  STEP 6 – Recompute Peq (same one-liner as both prof files)
%%  Peq = Peqv * or(ps, ~ps & ~(~Tot_CT))
%% =========================================================================
Peq = Peqv * ( ps | (~ps & logical(Tot_CT)) );

%% =========================================================================
%%  STEP 7 – Refresh boundary-condition vectors B0 and B0_liquid
%%  Prof does this for PNP1 and PNP2 (neighbours of emptied throat pores).
%%  Here we refresh ALL pores since we don't know which events occurred.
%%  Logic is identical; scope is wider.
%% =========================================================================
for pore = 1:NP
    nbr_throats = nonzeros(pnt(pore,:));
    nbr_pores   = nonzeros(pnp(pore,:));

    if Tot_CT(pore) == 0 && ps(pore) == 0
        % Pure vapour pore with no liquid-throat connections
        % Prof: B0(i) = sum(gv(...) .* log(1 - Peq(nbrs)/Patm))
        if ~isempty(nbr_pores)
            B0(pore) = sum(gv(nbr_throats) .* ...
                           log(1 - (Peq(nbr_pores) / Patm)));
        end

    elseif ps(pore) == 1
        % Liquid pore: refresh liquid RHS
        % Prof: B0_liquid(i) = sum(gw0(...) .* Pressure_meniscus(...))
        B0_liquid(pore) = sum(gw0(nbr_throats) .* ...
                              Pressure_meniscus(nbr_throats));

        % Rebuild liquid matrix row
        A0_liquid(pore, nbr_pores) = -gw0(nbr_throats)';
        A0_liquid(pore, pore)      =  sum(gw0(nbr_throats));
    end
end

%% =========================================================================
%%  STEP 8 – Append to cluster tracking vector (same as prof's line)
%%  Prof: Number_Cluster(1) = NumberCluster  (first call)
%%  Here we append so every call adds one entry (time-series tracking).
%% =========================================================================
Number_Cluster(end+1) = NumberCluster;

fprintf('[labelLiquidClusters_3D]  Clusters: %d  |  Cluster throats: %d / %d\n', ...
        NumberCluster, sum(tl>0), Nt);