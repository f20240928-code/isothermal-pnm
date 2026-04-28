%% =========================================================================
%%  solveLiquidPressure_3D.m
%%  Exact 3D replica of prof's  pressure-solution file
%%  (the file that starts:  "Peq = Peqv * or(ps, ~ps & ~(~Tot_CT))")
%%
%%  PURPOSE
%%  -------
%%  Given the current invasion state (ps, ts, gw0, Pressure_meniscus, …)
%%  this file:
%%    1. Decides which pores need equilibrium vapour pressure (Peq)
%%    2. Builds the liquid-phase sparse matrix  A0_liquid
%%    3. Populates the corresponding RHS  B0_liquid  and  B0 (vapour RHS)
%%    4. Computes per-pore total active cross-section  Total_Cross_A
%%    5. Counts clusters and initialises the saturation-slice tracker
%%
%%  CALL CONTEXT
%%  ------------
%%  Called inside the invasion loop AFTER solveVaporPressure_3D has run
%%  and after the invasion / throat-filling update step.
%%
%%  REQUIRES (workspace)
%%  --------------------
%%  network.static  : NP, pnp, pnt
%%  ps(NP,1)        : pore state  (1 = liquid-filled, 0 = vapour)
%%  Tot_CT(NP,1)    : liquid-throat coordination number  [from solveVaporPressure_3D]
%%  Peqv            : scalar equilibrium vapour pressure  [Pa]
%%  Patm            : atmospheric / reference pressure    [Pa]
%%  gv(NT,1)        : vapour conductance
%%  gw0(NT,1)       : liquid conductance
%%  Pressure_meniscus(NT,1): capillary pressure at menisci
%%  Atactive(NT,1)  : active cross-section areas
%%  tl(NT,…)        : throat-cluster label array
%%  Vt(Nt,1)        : throat volumes  (= network.static.Vt)
%%  m               : number of pore rows  (= network.static.m)
%%
%%  PRODUCES
%%  --------
%%  Peq(NP,1)        : equilibrium vapour pressure per pore
%%  A0_liquid        : NP×NP sparse liquid-pressure matrix
%%  B0_liquid(NP,1)  : RHS for liquid linear system
%%  B0(NP,1)         : RHS for vapour linear system (vapour-side entries)
%%  Volume_Liquid(Nt,1): current liquid volume per throat  (= Vt at entry)
%%  Total_Cross_A(NP,1): sum of active throat areas at each pore
%%  NumberCluster    : number of distinct liquid clusters
%%  Number_Cluster(1): same scalar, stored in tracking vector
%%  q                : saturation-profile save counter (initialised 99)
%%  Sslice(101,m-1)  : saturation profiles  (row 1 = all-liquid baseline)
%% =========================================================================

%% --- Unpack static network -----------------------------------------------
NP  = network.static.NP;
pnp = network.static.pnp;   % NP × 6
pnt = network.static.pnt;   % NP × 6
Vt  = network.static.Vt;    % Nt × 1
m   = network.static.m;

%% =========================================================================
%%  STEP 1 – Equilibrium vapour pressure  Peq
%%  2D:  Peq = Peqv * or(ps,  ~ps & ~(~Tot_CT))
%%
%%  Logic (unchanged from 2D):
%%    A pore gets Peq = Peqv  when it is either:
%%      (a) liquid-filled  [ps == 1], OR
%%      (b) vapour-filled but connected to at least one liquid throat
%%          [ps == 0  AND  Tot_CT > 0  →  ~(~Tot_CT) == logical(Tot_CT)]
%%
%%  In 3D, Tot_CT already sums over 6 neighbours (done in solveVaporPressure_3D),
%%  so the formula is IDENTICAL.
%% =========================================================================
Peq = Peqv * ( ps  |  (~ps & logical(Tot_CT)) );
%  Note: or(A,B) ≡ A|B in MATLAB; ~(~x) ≡ logical(x) for numeric x.

%% =========================================================================
%%  STEP 2 – Liquid-pressure sparse matrix  A0_liquid  and  RHS  B0_liquid
%%  2D:  spalloc(NP, NP, 5*NP)  → coordination 4+1
%%  3D:  spalloc(NP, NP, 7*NP)  → coordination 6+1
%% =========================================================================
Volume_Liquid = Vt;                          % initial liquid volume = full pore volume
A0_liquid     = spalloc(NP, NP, 7*NP);      % 3D: 7 = max coordination (6) + diagonal
B0_liquid     = zeros(NP, 1);

for i = 1:NP
    nbr_throats = nonzeros(pnt(i,:));        % up to 6 throat indices
    nbr_pores   = nonzeros(pnp(i,:));        % up to 6 neighbour pore indices

    % --- Off-diagonal: negative liquid conductance to each neighbour -----
    % 2D: A0_liquid(i, nonzeros(pnp(i,:))) = -gw0(nonzeros(pnt(i,:)))'
    A0_liquid(i, nbr_pores) = -gw0(nbr_throats)';

    % --- Diagonal: sum of outgoing liquid conductances -------------------
    % 2D: A0_liquid(i,i) = sum(gw0(nonzeros(pnt(i,:))))
    A0_liquid(i, i) = sum(gw0(nbr_throats));

    % --- RHS: liquid-filled pores ----------------------------------------
    % 2D:  if ps(i)==1
    %        B0_liquid(i) = sum(gw0(...) .* Pressure_meniscus(...))
    %      end
    if ps(i) == 1
        B0_liquid(i) = sum( gw0(nbr_throats) .* Pressure_meniscus(nbr_throats) );
    end

    % --- RHS: pure vapour pores (no Peq) ---------------------------------
    % 2D:  if ps(i)==0 & Peq(i)==0
    %        B0(i) = sum(gv(...) .* log(1 - Peq(nbr)/Patm))
    %      end
    % Note: B0 (vapour RHS) is updated here for vapour pores that have
    % no equilibrium pressure (isolated from liquid).
    if ps(i) == 0 && Peq(i) == 0
        nbr_pores_nz = nonzeros(pnp(i,:));   % need these for Peq indexing
        B0(i) = sum( gv(nbr_throats) .* log(1 - (Peq(nbr_pores_nz) / Patm)) );
    end
end

%% =========================================================================
%%  STEP 3 – Total active cross-sectional area per pore
%%  2D:  Total_Cross_A(i) = sum(Atactive(nonzeros(pnt(i,:))))
%%  3D:  identical formula; pnt now has 6 columns
%% =========================================================================
Total_Cross_A = zeros(NP, 1);
for i = 1:NP
    Total_Cross_A(i) = sum( Atactive(nonzeros(pnt(i,:))) );
end

%% =========================================================================
%%  STEP 4 – Cluster counting  and  saturation-slice tracker
%%  2D:  NumberCluster = max(max(tl))  → identical in 3D
%% =========================================================================
NumberCluster    = max(tl(:));            % max() over full array (works for any dim)
Number_Cluster(1) = NumberCluster;

% Saturation slice tracker
% 2D:  Sslice = zeros(101, m-1);  Sslice(1,:) = ones(1, m-1)
% 3D:  m-1 inter-layer gaps in the Y-direction; structure unchanged
q      = 99;                             % save every 1 % of total saturation
Sslice = zeros(101, m-1);
Sslice(1,:) = ones(1, m-1);             % row 1 = fully-liquid baseline

fprintf('[solveLiquidPressure_3D]  A0_liquid built: nnz=%d | Clusters=%d\n', ...
        nnz(A0_liquid), NumberCluster);