%% =========================================================================
%%  3D POROUS NETWORK GENERATOR
%%  Exact 3D physics replica of the 2D reference code
%%
%%  2D -> 3D mapping:
%%    2D: m x n  grid (rows x cols),  4 neighbours per pore (U/D/L/R)
%%    3D: m x n x p grid,             6 neighbours per pore (+X/-X/+Y/-Y/+Z/-Z)
%%
%%  Throat types (same taxonomy as 2D, extended):
%%    X-throats : connect (i,j,k)-(i,j+1,k)   [was "horizontal" in 2D]
%%    Y-throats : connect (i,j,k)-(i+1,j,k)   [was "vertical"   in 2D]
%%    Z-throats : connect (i,j,k)-(i,j,k+1)   [NEW third axis   in 3D]
%%    BL throats: boundary-layer connections   [same logic as 2D, now a slab]
%%
%%  Storage convention requested by user:
%%    network.static.*  – geometry / topology (never changes during simulation)
%%    network.dynamic.* – fields needed for video / time-stepping
%%
%%  Prerequisites – the following parameters must be defined BEFORE running:
%%    m, n, p      – interior grid dimensions (rows, cols, layers)
%%    mbl          – number of boundary-layer planes above the top face
%%    L4           – lattice spacing (interior)
%%    Lblh, Lblv   – boundary-layer horizontal / vertical spacings
%%    Ld           – depth used for boundary-layer throat areas
%%    rmean        – mean throat radius (interior)
%%  All other derived quantities are computed here.
%% =========================================================================

%% ------------------------------------------------------------------
%% 0.  DERIVED DIMENSIONS
%% ------------------------------------------------------------------
Np   = m * n * p;                % number of interior pores
Nsurf = n * p;                   % pores on one face (top XZ-plane)
NP   = Np + mbl * Nsurf;         % total pores  (interior + boundary layer)

% Interior throat counts
Ntx  = m * (n-1) * p;           % X-direction throats
Nty  = (m-1) * n * p;           % Y-direction throats (was "vertical" in 2D)
Ntz  = m * n * (p-1);           % Z-direction throats (new axis)
Nt   = Ntx + Nty + Ntz;         % total interior throats

% Boundary-layer throat counts  (same logic as 2D, face is n x p)
NtBLv = Nsurf * mbl;            % vertical BL connections (face->BL, BL->BL)
NtBLh = (n-1) * p * (mbl+1) + n * (p-1) * (mbl+1); % horizontal BL connections
NT    = Nt + NtBLv + NtBLh;     % grand total of throats

fprintf('3D Network: %d x %d x %d  |  Pores: %d  |  Throats: %d\n', m,n,p,NP,NT);

%% ------------------------------------------------------------------
%% 1.  PORE POSITIONS  (3-D Cartesian)
%% ------------------------------------------------------------------
pos = zeros(NP, 3);             % columns: X, Y, Z

% --- interior pores ---
for i = 1:m          % Y-index (rows, primary flow direction)
    for j = 1:n      % X-index (columns)
        for k = 1:p  % Z-index (depth)
            P = (i-1)*n*p + (j-1)*p + k;
            pos(P,1) = (j-1) * L4;   % X
            pos(P,2) = (i-1) * L4;   % Y
            pos(P,3) = (k-1) * L4;   % Z
        end
    end
end

% --- boundary layer pores  (slab above top Y-face) ---
for i = 1:mbl
    for j = 1:n
        for k = 1:p
            P = Np + (i-1)*Nsurf + (j-1)*p + k;
            pos(P,1) = (j-1) * Lblh;
            pos(P,2) = (m-1)*L4 + i * Lblv;   % above the top face
            pos(P,3) = (k-1) * Lblh;
        end
    end
end

%% ------------------------------------------------------------------
%% 2.  PNP – 6 NEIGHBOURS PER PORE
%%     Column order: [+Y , -X , -Y , +X , +Z , -Z]
%%                     up  left down right front back
%%     Zero = no neighbour (wall / edge)
%% ------------------------------------------------------------------
pnp = zeros(NP, 6);

% --- interior pores ---
for i = 1:m
    for j = 1:n
        for k = 1:p
            P = (i-1)*n*p + (j-1)*p + k;

            % +Y (up, i+1)
            if i ~= m
                pnp(P,1) = P + n*p;
            end
            % -X (left, j-1)
            if j ~= 1
                pnp(P,2) = P - p;
            end
            % -Y (down, i-1)
            if i ~= 1
                pnp(P,3) = P - n*p;
            end
            % +X (right, j+1)
            if j ~= n
                pnp(P,4) = P + p;
            end
            % +Z (front, k+1)
            if k ~= p
                pnp(P,5) = P + 1;
            end
            % -Z (back, k-1)
            if k ~= 1
                pnp(P,6) = P - 1;
            end
        end
    end
end

% --- boundary layer pores  (only Y and XZ lateral neighbours) ---
for i = 1:mbl
    for j = 1:n
        for k = 1:p
            P = Np + (i-1)*Nsurf + (j-1)*p + k;

            % +Y (up in BL)
            if i ~= mbl
                pnp(P,1) = P + Nsurf;
            end
            % -X
            if j ~= 1
                pnp(P,2) = P - p;
            end
            % -Y (down: either BL layer below or top face of interior)
            if i == 1
                % connect back to top interior face
                pnp(P,3) = (m-1)*n*p + (j-1)*p + k;   % top interior pore
            else
                pnp(P,3) = P - Nsurf;
            end
            % +X
            if j ~= n
                pnp(P,4) = P + p;
            end
            % +Z
            if k ~= p
                pnp(P,5) = P + 1;
            end
            % -Z
            if k ~= 1
                pnp(P,6) = P - 1;
            end
        end
    end
end

pnp = -sort(-pnp, 2);           % zeros pushed to the right (same as 2D)

%% ------------------------------------------------------------------
%% 3.  TNP – TWO NEIGHBOURING PORES PER THROAT
%% ------------------------------------------------------------------
tnp = zeros(NT, 2);

% --- Y-direction throats (rows, was "vertical" in 2D) ---
idx = 0;
for i = 1:m-1
    for j = 1:n
        for k = 1:p
            idx = idx + 1;
            P = (i-1)*n*p + (j-1)*p + k;
            tnp(idx, 1) = P;
            tnp(idx, 2) = P + n*p;
        end
    end
end
% quick check: idx should equal Nty
assert(idx == Nty, 'Y-throat count mismatch');

% --- X-direction throats (columns, was "horizontal" in 2D) ---
for i = 1:m
    for j = 1:n-1
        for k = 1:p
            idx = idx + 1;
            P = (i-1)*n*p + (j-1)*p + k;
            tnp(idx, 1) = P;
            tnp(idx, 2) = P + p;
        end
    end
end
assert(idx == Nty + Ntx, 'X-throat count mismatch');

% --- Z-direction throats (depth, new axis) ---
for i = 1:m
    for j = 1:n
        for k = 1:p-1
            idx = idx + 1;
            P = (i-1)*n*p + (j-1)*p + k;
            tnp(idx, 1) = P;
            tnp(idx, 2) = P + 1;
        end
    end
end
assert(idx == Nt, 'Interior throat count mismatch');

% --- Boundary-layer vertical connections (face -> BL -> BL) ---
for i = 1:mbl
    for j = 1:n
        for k = 1:p
            idx = idx + 1;
            P_bl = Np + (i-1)*Nsurf + (j-1)*p + k;
            if i == 1
                P_below = (m-1)*n*p + (j-1)*p + k;  % top interior face
            else
                P_below = P_bl - Nsurf;
            end
            tnp(idx, 1) = P_below;
            tnp(idx, 2) = P_bl;
        end
    end
end
assert(idx == Nt + NtBLv, 'BL vertical throat count mismatch');

% --- Boundary-layer horizontal connections (X and Z directions) ---
for i = 0:mbl        % i=0 is the top interior face itself
    for j = 1:n
        for k = 1:p-1   % Z-direction
            idx = idx + 1;
            if i == 0
                P = (m-1)*n*p + (j-1)*p + k;
            else
                P = Np + (i-1)*Nsurf + (j-1)*p + k;
            end
            tnp(idx, 1) = P;
            tnp(idx, 2) = P + 1;
        end
    end
    for j = 1:n-1       % X-direction
        for k = 1:p
            idx = idx + 1;
            if i == 0
                P = (m-1)*n*p + (j-1)*p + k;
            else
                P = Np + (i-1)*Nsurf + (j-1)*p + k;
            end
            tnp(idx, 1) = P;
            tnp(idx, 2) = P + p;
        end
    end
end
assert(idx == NT, sprintf('Total throat count mismatch: got %d expected %d', idx, NT));

%% ------------------------------------------------------------------
%% 4.  PNT – THROAT INDEX AT EACH PORE NEIGHBOUR SLOT
%% ------------------------------------------------------------------
pnt = zeros(NP, 6);
for j = 1:NT
    k = tnp(j,1);
    r = tnp(j,2);
    slot_k = find(pnp(k,:) == r, 1);
    slot_r = find(pnp(r,:) == k, 1);
    if ~isempty(slot_k), pnt(k, slot_k) = j; end
    if ~isempty(slot_r), pnt(r, slot_r) = j; end
end

%% ------------------------------------------------------------------
%% 5.  TNT / TNP2 – THROAT NEIGHBOUR THROATS (same logic as 2D)
%% ------------------------------------------------------------------
tnt  = zeros(Nt, 10);   % up to 10 neighbour throats in 3D (2*5 max)
tnp2 = zeros(Nt, 10);

for k = 1:Nt
    g1   = nonzeros(pnt(tnp(k,1), :));
    g1   = g1(g1 ~= k)';
    Nh1  = length(g1);
    h1   = tnp(k,1) * ones(1, Nh1);

    g2   = nonzeros(pnt(tnp(k,2), :));
    g2   = g2(g2 ~= k)';
    Nh2  = length(g2);
    h2   = tnp(k,2) * ones(1, Nh2);

    tnt (k, 1:Nh1+Nh2) = [g1  g2];
    tnp2(k, 1:Nh1+Nh2) = [h1  h2];
end

%% ------------------------------------------------------------------
%% 6.  RADII – 3D deterministic pattern with enhanced cross-planes
%%     Mirrors the 2D sin/cos pattern; adds Z-variation for the new axis.
%%     "Cross" pattern: middle X-column, middle Y-row, middle Z-layer enlarged.
%% ------------------------------------------------------------------
rt = zeros(Nt, 1);

% Y-direction throats  (analogue of 2D "vertical")
for i = 1:m-1
    for j = 1:n
        for k = 1:p
            t = (i-1)*n*p + (j-1)*p + k;
            variation = 1 + 0.05 * sin(i*pi/m) * cos(j*pi/n) * cos(k*pi/p);
            rt(t) = rmean * variation;
        end
    end
end

% X-direction throats  (analogue of 2D "horizontal")
base = Nty;
for i = 1:m
    for j = 1:n-1
        for k = 1:p
            t = base + (i-1)*(n-1)*p + (j-1)*p + k;
            variation = 1 + 0.05 * cos(i*pi/m) * sin(j*pi/n) * cos(k*pi/p);
            rt(t) = rmean * variation;
        end
    end
end

% Z-direction throats  (new)
base = Nty + Ntx;
for i = 1:m
    for j = 1:n
        for k = 1:p-1
            t = base + (i-1)*n*(p-1) + (j-1)*(p-1) + k;
            variation = 1 + 0.05 * cos(i*pi/m) * cos(j*pi/n) * sin(k*pi/p);
            rt(t) = rmean * variation;
        end
    end
end

% Enlarge middle X-column (analogue of 2D middle column)
mid_j = round(n/2);
for i = 1:m-1
    for k = 1:p
        t = (i-1)*n*p + (mid_j-1)*p + k;
        rt(t) = rt(t) * 2.5;
    end
end

% Enlarge middle Y-row (analogue of 2D middle row)
mid_i = round(m/2);
for j = 1:n-1
    for k = 1:p
        t = Nty + (mid_i-1)*(n-1)*p + (j-1)*p + k;
        rt(t) = rt(t) * 2.5;
    end
end

% Enlarge middle Z-layer (new cross-plane, no 2D analogue)
mid_k = round(p/2);
for i = 1:m
    for j = 1:n
        t = Nty + Ntx + (i-1)*n*(p-1) + (j-1)*(p-1) + mid_k;
        if mid_k <= p-1
            rt(t) = rt(t) * 2.5;
        end
    end
end

save('radius_3D', 'rt');

%% ------------------------------------------------------------------
%% 7.  THROAT AREAS, LENGTHS, VOLUMES
%% ------------------------------------------------------------------
At = zeros(NT, 1);

% Interior throats: circular cross-section
At(1:Nt) = pi * rt .* rt;

% Boundary-layer throats (rectangular, same formula as 2D extended to 3D)
At(Nt+1 : Nt+NtBLv)            = Ld * Lblh * ones(NtBLv, 1);
At(Nt+NtBLv+1 : NT)            = Ld * Lblv * ones(NtBLh, 1);

% Lengths
Lt = [L4*ones(Nt,1); Lblv*ones(NtBLv,1); Lblh*ones(NtBLh,1)];

% Volumes (interior throats only, same as 2D)
Vt = At(1:Nt) .* Lt(1:Nt);

% Slice volumes  (Y-direction slices, for flux diagnostics – same as 2D Vtslice)
Nty_per_slice = n * p;               % Y-throats per inter-layer gap
Ntx_per_slice = (n-1) * p;          % X-throats per layer
Ntz_per_slice = n * (p-1);          % Z-throats per layer

Vtslice = zeros(1, m-1);
for sl = 1:m-1
    yT  = Vt( (sl-1)*Nty_per_slice + 1 : sl*Nty_per_slice );
    xT  = Vt( Nty + (sl-1)*Ntx_per_slice + 1 : Nty + sl*Ntx_per_slice );
    zT  = Vt( Nty+Ntx + (sl-1)*Ntz_per_slice + 1 : Nty+Ntx + sl*Ntz_per_slice );
    Vtslice(sl) = sum(yT) + sum(xT) + sum(zT);
end

TotalVolume_Pore_Network = sum(Vt);
fprintf('Total interior pore-network volume = %.6e\n', TotalVolume_Pore_Network);

%% ------------------------------------------------------------------
%% 8.  PACK INTO network.static  AND  network.dynamic
%% ------------------------------------------------------------------

% ---- STATIC: topology + geometry (never mutated during simulation) ----
network.static.NP       = NP;
network.static.Np       = Np;
network.static.NT       = NT;
network.static.Nt       = Nt;
network.static.Ntx      = Ntx;
network.static.Nty      = Nty;
network.static.Ntz      = Ntz;
network.static.NtBLv    = NtBLv;
network.static.NtBLh    = NtBLh;
network.static.m        = m;
network.static.n        = n;
network.static.p        = p;
network.static.mbl      = mbl;
network.static.Nsurf    = Nsurf;
network.static.L4       = L4;
network.static.Lblh     = Lblh;
network.static.Lblv     = Lblv;
network.static.Ld       = Ld;

% connectivity
network.static.pos      = pos;    % NP x 3  Cartesian coordinates
network.static.pnp      = pnp;    % NP x 6  neighbour pore indices
network.static.tnp      = tnp;    % NT x 2  pore pair per throat
network.static.pnt      = pnt;    % NP x 6  throat index at each neighbour slot
network.static.tnt      = tnt;    % Nt x 10 neighbouring-throat indices
network.static.tnp2     = tnp2;   % Nt x 10 shared-pore indices for tnt

% geometry
network.static.rt       = rt;     % Nt x 1  interior throat radii
network.static.At       = At;     % NT x 1  throat cross-sectional areas
network.static.Lt       = Lt;     % NT x 1  throat lengths
network.static.Vt       = Vt;     % Nt x 1  throat volumes
network.static.Vtslice  = Vtslice;% 1 x m-1 slice volumes (diagnostics)
network.static.TotalVolume = TotalVolume_Pore_Network;

% ---- DYNAMIC: fields written / read every time-step for video ----
% Initialised to zero here; the flow/transport solver updates these.
network.dynamic.Pressure     = zeros(NP, 1);   % pore pressures
network.dynamic.Saturation   = zeros(NP, 1);   % fluid saturation per pore
network.dynamic.FlowRate     = zeros(NT, 1);   % volumetric flow per throat
network.dynamic.Velocity     = zeros(NT, 1);   % mean velocity per throat
network.dynamic.Occupancy    = zeros(NT, 1);   % 0=gas 1=liquid (invasion flag)
network.dynamic.Concentration= zeros(NP, 1);   % solute/tracer concentration
network.dynamic.Time         = 0;              % current simulation time
network.dynamic.Step         = 0;              % current time-step counter

% ---- SAVE ----
save('network_3D.mat', 'network', '-v7.3');
fprintf('Saved: network_3D.mat  (network.static + network.dynamic)\n');

%% ------------------------------------------------------------------
%% 9.  QUICK SANITY PLOT  (3-D wire-frame of interior throats)
%% ------------------------------------------------------------------
figure(1); clf; hold on;
% Y-direction throats – blue
for t = 1:Nty
    P1 = tnp(t,1); P2 = tnp(t,2);
    line([pos(P1,1) pos(P2,1)],[pos(P1,2) pos(P2,2)],[pos(P1,3) pos(P2,3)],...
        'Color','b','LineWidth',0.5);
end
% X-direction throats – black
for t = Nty+1:Nty+Ntx
    P1 = tnp(t,1); P2 = tnp(t,2);
    line([pos(P1,1) pos(P2,1)],[pos(P1,2) pos(P2,2)],[pos(P1,3) pos(P2,3)],...
        'Color','k','LineWidth',0.5);
end
% Z-direction throats – green
for t = Nty+Ntx+1:Nt
    P1 = tnp(t,1); P2 = tnp(t,2);
    line([pos(P1,1) pos(P2,1)],[pos(P1,2) pos(P2,2)],[pos(P1,3) pos(P2,3)],...
        'Color',[0 0.6 0],'LineWidth',0.5);
end
% Boundary layer – red
for t = Nt+1:NT
    P1 = tnp(t,1); P2 = tnp(t,2);
    line([pos(P1,1) pos(P2,1)],[pos(P1,2) pos(P2,2)],[pos(P1,3) pos(P2,3)],...
        'Color','r','LineWidth',0.8);
end
% Pore markers
plot3(pos(1:Np,1),    pos(1:Np,2),    pos(1:Np,3),    'ko','MarkerFaceColor','k','MarkerSize',2);
plot3(pos(Np+1:end,1),pos(Np+1:end,2),pos(Np+1:end,3),'ro','MarkerFaceColor','w','MarkerSize',2);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Porous Network (blue=Y, black=X, green=Z throats, red=BL)');
axis equal; grid on; view(30,25);