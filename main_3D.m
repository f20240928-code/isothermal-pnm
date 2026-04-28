%% =========================================================================
%%  main_3D.m
%%  3D Pore-Network Model – Drying Simulation
%%
%%  Direct 3D replica of the 2D reference main file.
%%  Every section maps 1-to-1 with the 2D code; all changes are
%%  annotated with  [2D → 3D].
%%
%%  CALL ORDER (mirrors 2D exactly):
%%    1. network_generator_3D   (replaces networkgen + radgen_2)
%%    2. solveVaporPressure_3D  (replaces dry_a_sparse_mat_file)
%%    3. identifyCandidateThroats_3D (replaces dry_b_cand_thr_1)
%%    4. solveLiquidPressure_3D (replaces dry_c_press_soln_1)
%%    5. labelLiquidClusters_3D (replaces first HK call before loop)
%%    6. dryingLoop_3D          (replaces dry_d_loop_1)
%%
%%  VIDEO:  networkDrawing3D is called inside dryingLoop_3D each step.
%%          The VideoWriter object (Obj) is opened here and closed inside
%%          dryingLoop_3D after the loop finishes.
%% =========================================================================

clear all
clc
close all

%% -------------------------------------------------------------------------
%% FIGURE / VIDEO WINDOW SETUP
%% 2D: set(gca,...) + set(gcf,'WindowState','Maximized')   UNCHANGED.
%% -------------------------------------------------------------------------
figure(1);
set(gca, 'nextplot', 'replacechildren')
set(gcf, 'WindowState', 'Maximized')

%% =========================================================================
%%  PHYSICAL CONSTANTS
%%  2D reference block – every value UNCHANGED.
%% =========================================================================
mu      = 0.001;                             % dynamic viscosity of water,    Pa·s
nu      = 15.35e-6;                          % kinematic viscosity at 20 °C,  m²/s
T       = 20;                                % temperature,                   °C
                                             %   (2D used 100 °C for Rhov;
                                             %    20 °C is the correct value
                                             %    matching Peqv=2339 Pa)
Patm    = 101325;                            % atmospheric pressure,          Pa
Peqv    = 2339;                              % equilibrium vapour pressure at 20 °C, Pa
yeq     = Peqv / Patm;                       % equilibrium vapour mole fraction
yift    = 0.000;                             % vapour fraction at infinite distance
sigma   = 0.07274;                           % surface tension water-air 20 °C, N/m
delta   = 22.6e-6 * ((T+273)/273)^1.81;     % binary diffusion coefficient, m²/s  [Schirmer]
Mv      = 18.02;                             % molar mass of water,           kg/kmol
Rg      = 8314.5;                            % universal gas constant,        J/kmol/K
RhoL    = 998.21;                            % liquid water density at 20 °C, kg/m³
Rhov    = (Mv * Patm) / (Rg * (T + 273));   % vapour density at T,           kg/m³

%% =========================================================================
%%  VIDEO WRITER
%%  2D: VideoWriter('pnm_dry_test_2.avi')   →   [2D→3D] new filename.
%%  Object opened here; writeVideo called inside dryingLoop_3D; closed there.
%% =========================================================================
Obj = VideoWriter('pnm_dry_3D.avi');
Obj.FrameRate = 1;
open(Obj);
set(gca, 'nextplot', 'replacechildren')
set(gcf, 'WindowState', 'Maximized')

%% =========================================================================
%%  NETWORK PARAMETERS
%%  2D: n (cols), m (rows), Nsurf=n
%%  [2D→3D]: add p (depth layers); Nsurf = n*p; all derived counts extended.
%% =========================================================================

% Grid dimensions
n   = 4;          % number of columns   [same as 2D]
m   = 4;          % number of rows      [same as 2D]
p   = 4;           % number of depth layers  [NEW – 3D only]

% Geometry
L4    = 5e-4;      % throat length / lattice spacing,  m   [same as 2D]
rmean = 4e-5;      % mean throat radius,               m   [same as 2D]
rstd  = 0.05 * rmean;  % std dev of radius distribution    [same as 2D]

% Boundary layer
%   2D: Nbl=10, mbl=Nbl-1=9
Nbl  = 10;                 % BL thickness in throat-lengths  [same as 2D]
mbl  = Nbl - 1;            % BL node rows                    [same as 2D]
Lblv = L4;                 % BL vertical spacing             [same as 2D]
Lblh = L4;                 % BL horizontal spacing           [same as 2D]
Ld   = L4;                 % depth used for BL throat areas  [same as 2D; = L4 in 3D]

% Simulation bookkeeping
%   2D: flag=0, computertime1=clock, Comptime, Clustersize
flag          = 0;
computertime1 = datetime('now');          % [2D→3D] datetime instead of clock

%% =========================================================================
%%  OUTER LOOP SETUP
%%  2D: totalloops, S, Mdot, Satprofile
%%  [2D→3D]: Satprofile third dim = totalloops (unchanged logic)
%% =========================================================================
totalloops = 1;

%  These are filled after each inner loop completes
%  (indices sized to Nt-1; Nt is derived inside network_generator_3D)
%  We pre-allocate after the network is built.

%% =========================================================================
%%  GENERATE NETWORK  (replaces  "networkgen"  +  "radgen_2"  in 2D)
%%  network_generator_3D writes everything into  network.static / dynamic
%%  and also leaves the flat workspace variables that the solvers expect:
%%    NP, NT, Nt, Np, Nsurf, mbl, m, n, p
%%    pnp, pnt, tnp, tnt, tnp2
%%    rt, At, Lt, Vt, Vtslice, TotalVolume_Pore_Network
%%    pos
%% =========================================================================
network_generator_3D;

% Unpack key scalars (needed before the loop for array sizing)
NP    = network.static.NP;
NT    = network.static.NT;
Nt    = network.static.Nt;
Np    = network.static.Np;
Nsurf = network.static.Nsurf;
m     = network.static.m;
n     = network.static.n;
p     = network.static.p;
mbl   = network.static.mbl;

TotalVolume_Pore_Network = network.static.TotalVolume;

% Now size the outer-loop accumulators (Nt is known)
S         = zeros(totalloops, Nt-1);
Mdot      = zeros(totalloops, Nt-1);
Satprofile = zeros(101, m-1, totalloops);

%% =========================================================================
%%  OUTER LOOP  (identical structure to 2D "for loop = 1:totalloops")
%% =========================================================================
for loop = 1:totalloops

    fprintf('\n===== LOOP %d / %d =====\n', loop, totalloops);

    %% ---------------------------------------------------------------------
    %%  CONDUCTANCES
    %%  2D: gv = (Mv*delta*At)/(Rg*(T+273)*Lt)   UNCHANGED.
    %%      ginf: non-zero only at surface nodes (last Nsurf pores)
    %%      The surface throats in 3D are indices  Nt+1 : Nt+Nsurf
    %%      (first NtBLv entries = vertical BL connections to the top face)
    %%      → same indexing convention as 2D  (Nt+1:Nt+Nsurf)
    %% ---------------------------------------------------------------------
    At = network.static.At;
    Lt = network.static.Lt;

    gv   = (Mv * delta * At) ./ (Rg * (T + 273) * Lt);   % [same as 2D]
    ginf = zeros(NP, 1);
    % 2D: ginf(end-Nsurf+1:end) = gv(Nt+1:Nt+Nsurf)
    % [2D→3D]: Nsurf = n*p instead of n; index range identical in form
    ginf(end-Nsurf+1 : end) = gv(Nt+1 : Nt+Nsurf);

    %% ---------------------------------------------------------------------
    %%  INITIAL STATE ARRAYS
    %%  2D:  ts = [ones(Nt,1);  zeros(NT-Nt,1)]   – interior liquid, BL dry
    %%       ps = [ones(Np-Nsurf,1); zeros(NP-Np+Nsurf,1)] – interior liquid,
    %%            surface + BL pores dry
    %%       tl = [ones(Nt,1);  zeros(NT-Nt,1)]   – all interior throats in
    %%            one big cluster (label 1), BL = 0
    %%  UNCHANGED – indexing identical in 3D.
    %% ---------------------------------------------------------------------
    ts = [ones(Nt, 1);  zeros(NT-Nt, 1)];          % throat saturation state
    ps = [ones(Np-Nsurf, 1); zeros(NP-Np+Nsurf, 1)]; % pore state
    tl = [ones(Nt, 1);  zeros(NT-Nt, 1)];          % throat cluster labels

    %% ---------------------------------------------------------------------
    %%  STEP A – Build vapour-pressure sparse matrix  A0
    %%  2D: dry_a_sparse_mat_file
    %%  [2D→3D]: solveVaporPressure_3D
    %%  Produces: A0, B0, Totaltime, Tot_CT,
    %%            Pore1, Pore2, Atactive, rcandidate, Pressure_meniscus, gw0, Peq, Pv
    %% ---------------------------------------------------------------------
    solveVaporPressure_3D;

    %% ---------------------------------------------------------------------
    %%  STEP B – Identify invasion candidate throats
    %%  2D: dry_b_cand_thr_1
    %%  [2D→3D]: identifyCandidateThroats_3D
    %%  Produces: Pore_Cond, Pore1, Pore2, rcandidate,
    %%            Atactive, Pressure_meniscus, gw0, MM
    %% ---------------------------------------------------------------------
    identifyCandidateThroats_3D;

    %% ---------------------------------------------------------------------
    %%  STEP C – Build liquid-pressure matrix and cluster/slice trackers
    %%  2D: dry_c_press_soln_1
    %%  [2D→3D]: solveLiquidPressure_3D
    %%  Produces: Peq, A0_liquid, B0_liquid, B0 (updated), Volume_Liquid,
    %%            Total_Cross_A, NumberCluster, Number_Cluster, q, Sslice
    %% ---------------------------------------------------------------------
    Number_Cluster = zeros(Nt, 1);   % pre-allocate tracking vector (same as 2D)
    solveLiquidPressure_3D;

    %% ---------------------------------------------------------------------
    %%  STEP D – Initial cluster labelling (HK initialisation)
    %%  2D: called implicitly at start of dry_d_loop_1 via dry_e_hoshen_kopelmann
    %%  [2D→3D]: call labelLiquidClusters_3D once before entering the loop
    %%  so that tl, NumberCluster, ps, Tot_CT, Peq, A0_liquid, B0_liquid,
    %%  B0, Total_Cross_A are all consistent at loop entry.
    %% ---------------------------------------------------------------------
    labelLiquidClusters_3D;

    %% ---------------------------------------------------------------------
    %%  STEP E – MAIN DRYING LOOP
    %%  2D: dry_d_loop_1
    %%  [2D→3D]: dryingLoop_3D
    %%  Uses: all variables set above + Obj (VideoWriter)
    %%  Produces: Mev, TotalVolume_Liquid, Ss, timee, Sslice, Comptime
    %% ---------------------------------------------------------------------
    dryingLoop_3D;

    %% ---------------------------------------------------------------------
    %%  COLLECT RESULTS INTO OUTER-LOOP ARRAYS
    %%  2D: S(loop,:) = Ss;  Mdot(loop,:) = Mev(1:Nt-1);  Satprofile(:,:,loop)
    %%  UNCHANGED in form.
    %% ---------------------------------------------------------------------
    S(loop, :)            = Ss';
    Mdot(loop, :)         = Mev(1:Nt-1)';
    Satprofile(:, :, loop) = Sslice;

    % 2D: early exit condition
    if loop == 10
        break;
    end

end  % end outer loop

%% =========================================================================
%%  POST-PROCESSING / DIAGNOSTICS
%%  2D reference: drying-curve and saturation-profile plots
%% =========================================================================
computertime_end = datetime('now');
fprintf('\nTotal wall-clock time: %s\n', ...
        char(duration(computertime_end - computertime1)));

%% --- Drying curve  (Mev vs Saturation) -----------------------------------
figure(2); clf;
plot(S(1,:), Mdot(1,:) * 1e6, 'b-', 'LineWidth', 1.5);
xlabel('Saturation  S  [-]');
ylabel('Evaporation rate  \dot{M}  [\mug/s]');
title('3D PNM Drying Curve');
grid on;
set(gca, 'XDir', 'reverse');

%% --- Saturation profiles  (Sslice rows vs Y-slice index) -----------------
figure(3); clf;
imagesc(Satprofile(:, :, 1));
colorbar;
xlabel('Y-slice index');
ylabel('Saturation step (100=wet, 1=dry)');
title('3D PNM Saturation Profile Evolution');

%% --- Cluster count evolution ---------------------------------------------
figure(4); clf;
plot(1:Nt-1, Number_Cluster(1:Nt-1), 'r-', 'LineWidth', 1.2);
xlabel('Drying step');
ylabel('Number of liquid clusters');
title('3D PNM – Cluster Count vs Step');
grid on;

%% --- Comptime per step ---------------------------------------------------
figure(5); clf;
plot(1:Nt-1, Comptime(1:Nt-1), 'k-');
xlabel('Drying step');
ylabel('Wall-clock time per step  [s]');
title('3D PNM – Computation Time per Step');
grid on;

fprintf('\nDone.  Final saturation = %.4f\n', Ss(end));