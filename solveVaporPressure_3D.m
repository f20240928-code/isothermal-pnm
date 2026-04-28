%% =========================================================================
%%  solveVaporPressure_3D.m
%%  Exact 3D replica of prof's  dry_sparse_mat_file.m
%%
%%  PURPOSE
%%  -------
%%  Builds the sparse linear system  A·Pv = B  for the vapour-phase pores
%%  and initialises every working array that the invasion / wetting loop
%%  will need.
%%
%%  REQUIRES (from workspace / main file)
%%  --------------------------------------
%%  network.static  : NP, NT, Nt, pnp, pnt, pos, At, Lt, rt  (6-neighbour)
%%  gv(NT,1)        : vapour conductance of every throat
%%  ginf(NP,1)      : external vapour-inlet conductance at each pore
%%  ts(NT,1)        : throat state  (1 = liquid-filled, 0 = vapour)
%%
%%  PRODUCES (added to workspace for use by solveLiquidPressure_3D)
%%  ----------------------------------------------------------------
%%  Totaltime       : simulation clock, reset to 0 here
%%  Tot_CT(NP,1)    : active liquid-throat coordination number per pore
%%  A0              : NP×NP sparse vapour-pressure matrix
%%  Pv(NP,1)        : vapour pressure vector (initialised 0, solved in loop)
%%  Peq(NP,1)       : equilibrium vapour pressure (initialised 0)
%%  B0(NP,1)        : RHS vector for vapour system
%%  Pore1/Pore2     : pore indices for each throat (working arrays)
%%  Atactive(NT,1)  : active (liquid) cross-sectional area per throat
%%  rcandidate(NT,1): invasion-candidate radii
%%  Pressure_meniscus(NT,1): capillary pressure at each meniscus
%%  gw0(NT,1)       : liquid conductance (initialised 0, filled in loop)
%% =========================================================================

%% --- Unpack static network (keeps code identical in style to 2D) --------
NP   = network.static.NP;
NT   = network.static.NT;
pnp  = network.static.pnp;   % NP × 6   (was NP × 4 in 2D)
pnt  = network.static.pnt;   % NP × 6   (was NP × 4 in 2D)

%% --- Simulation clock ----------------------------------------------------
Totaltime = 0;

%% --- Tot_CT : sum of liquid-throat states at each pore ------------------
%  Exact replica of 2D line:
%    Tot_CT(i) = sum(ts(nonzeros(pnt(i,:))))
%  The only change: pnt now has 6 columns instead of 4.
Tot_CT = zeros(NP, 1);

%% --- Sparse matrix A0 for vapour pressures ------------------------------
%  2D used spalloc(NP, NP, 5*NP)  where 5 = coordination + 1.
%  3D coordination number = 6  →  allocate  7*NP  non-zeros.
A0 = spalloc(NP, NP, 7*NP);

for i = 1:NP
    nbr_throats = nonzeros(pnt(i,:));          % active throat indices (up to 6)
    nbr_pores   = nonzeros(pnp(i,:));          % neighbour pore indices

    % Off-diagonal: negative conductance to each neighbour
    % 2D: A0(i, nonzeros(pnp(i,:))) = -gv(nonzeros(pnt(i,:)))'
    A0(i, nbr_pores) = -gv(nbr_throats)';

    % Diagonal: sum of all outflow conductances + external vapour inlet
    % 2D: A0(i,i) = sum(gv(nonzeros(pnt(i,:)))) + ginf(i)
    A0(i, i) = sum(gv(nbr_throats)) + ginf(i);

    % Active liquid-throat coordination number (same formula, 6-slot pnt)
    % 2D: Tot_CT(i) = sum(ts(nonzeros(pnt(i,:))))
    Tot_CT(i) = sum(ts(nbr_throats));
end

%% --- Initialise all working vectors -------------------------------------
%  Exact mirrors of the 2D initialisations.

Pv   = zeros(NP, 1);    % vapour pressure at each pore
Peq  = zeros(NP, 1);    % equilibrium vapour pressure at each pore
B0   = zeros(NP, 1);    % RHS vector for vapour linear system

Pore1 = zeros(NT, 1);   % pore-1 index for each throat (filled during invasion)
Pore2 = zeros(NT, 1);   % pore-2 index for each throat

Atactive          = zeros(NT, 1);  % active (liquid) cross-section area
rcandidate        = zeros(NT, 1);  % invasion-candidate throat radii
Pressure_meniscus = zeros(NT, 1);  % capillary pressure at each meniscus
gw0               = zeros(NT, 1);  % liquid conductance (Hagen-Poiseuille)

fprintf('[solveVaporPressure_3D]  A0 built: %d pores, %d throats, nnz(A0)=%d\n', ...
        NP, NT, nnz(A0));