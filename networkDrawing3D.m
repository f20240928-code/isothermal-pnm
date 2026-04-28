%% =========================================================================
%%  networkDrawing3D.m
%%  3D replica of the prof's 2D network drawing / video-frame script.
%%
%%  2D DESIGN (prof's code)
%%  -----------------------
%%  subplot(1,2,1):  2D plan-view wire diagram of the network.
%%    - Grey rectangle = network bounding box
%%    - Throat colour:
%%        black (ts==1)  → fully liquid
%%        white (ts==0)  → fully vapour
%%        blue           → partially filled  (0 < ts < 1)
%%    - Line width proportional to throat radius (replaces "large" list)
%%  subplot(1,2,2):  live drying curve  Mev vs Ss
%%    - x-axis  Ss  (saturation, 0→1),  y-axis  Mev  [kg/s]
%%
%%  3D EXTENSION
%%  ------------
%%  The 3D network cannot be faithfully shown in a single 2D plan view.
%%  Instead we use THREE subplots that preserve EVERY visual convention
%%  of the prof exactly:
%%
%%    subplot(1,3,1)  – XY projection  (front view, col vs row)
%%                      Direct analogue of the prof's subplot(1,2,1).
%%                      Grey box, same colour/width rules, same title.
%%    subplot(1,3,2)  – XZ projection  (top view, col vs depth)
%%                      Same colour/width rules; complements the front view.
%%    subplot(1,3,3)  – Live drying curve  Mev vs Ss
%%                      Exact replica of the prof's subplot(1,2,2).
%%
%%  COLOUR CONVENTION  (identical to 2D)
%%    black  → ts == 1  (fully liquid)
%%    white  → ts == 0  (fully vapour / empty)
%%    blue   → 0 < ts < 1  (meniscus / partially filled)
%%
%%  LINE WIDTH CONVENTION  (replaces the prof's hard-coded "large" list)
%%    lw = lw_base + (lw_scale * rt(k)/rmean)
%%    Interior throats only (k ≤ Nt); BL throats use lw_base.
%%    This is a physics-based analogue of the prof's large-throat emphasis.
%%
%%  REQUIRES (workspace – all set by dryingLoop_3D at each step)
%%  -------------------------------------------------------------
%%  network.static : pos [NP×3], tnp [NT×2], rt [Nt×1]
%%                   NT, Nt, m, n, p, Nsurf, L4, Lblv, rmean
%%  ts(NT,1)       : throat saturation  (1=full liq, 0=vapour, else partial)
%%  Ss(Nt-1,1)     : saturation history  (Ss(s) set at current step s)
%%  Mev(NT,1)      : evaporation-rate history  (Mev(s) set at current step)
%%  s              : current drying-loop step index
%%
%%  CALLED BY
%%  ---------
%%  dryingLoop_3D.m  (inside the for s = 1:Nt-1 loop, just before writeVideo)
%% =========================================================================

%% --- Unpack static fields -----------------------------------------------
pos    = network.static.pos;     % NP × 3  [X Y Z]
tnp    = network.static.tnp;     % NT × 2
rt     = network.static.rt;      % Nt × 1  (interior throats only)
NT     = network.static.NT;
Nt     = network.static.Nt;
m      = network.static.m;
n      = network.static.n;
p      = network.static.p;
L4     = network.static.L4;
Lblv   = network.static.Lblv;
rmean_ = mean(rt);               % local mean for line-width scaling

%% --- Classify interior throats (indices 1:Nt) ---------------------------
%  Same logic as 2D:
%    blacks = ts == 1  (fully liquid)
%    whites = ts == 0  (fully vapour)
%    blues  = everything else (partial fill / meniscus)
%  Only interior throats (1:Nt) are drawn in the phase panels, identical
%  to the prof who loops over 1:Nt.  BL throats are omitted (same as 2D).

blacks = find(ts(1:Nt) == 1);
whites = find(ts(1:Nt) == 0);
blues  = setdiff(1:Nt, union(blacks, whites));

%% --- Line-width helper --------------------------------------------------
%  Prof used a hard-coded "large" index list and set Linewidth=5 for those.
%  Here we scale continuously with radius (same physical intent):
%    lw = 0.8  (thin)  for rt ≈ rmean
%    lw = 2.5  (thick) for enlarged throats  (rt up to 2.5*rmean)
lw_base  = 0.5;
lw_scale = 2.0;    % extra width per (rt/rmean) above baseline

lw_for = @(k) lw_base + lw_scale * (rt(k) / rmean_);   % k = interior index

%% =========================================================================
%%  SUBPLOT 1  –  XY PROJECTION  (front view, replicates prof's subplot(1,2,1))
%%
%%  Axis: X = pos(:,1),  Y = pos(:,2)
%%  Grey bounding box identical to prof:
%%    Rectangle([-L4/2, -L4/2,  width,  height])
%%    width  = L4*n   (same as 2D where n = number of columns)
%%    height = L4*(m-0.5)
%% =========================================================================
subplot(1, 3, 1);

% --- Grey background box (exact replica of prof's rectangle) -------------
rectangle('Position', [-L4/2, -L4/2, L4*n, L4*(m-0.5)], ...
          'Curvature', [0 0], ...
          'FaceColor', [0.6 0.6 0.6], ...
          'EdgeColor', 'none');
hold on;

% --- Blue throats (partial fill) -----------------------------------------
%  2D: if ~isempty(blues) → for z=1:length(blues) → line(...,'color','b')
if ~isempty(blues)
    for z = 1:length(blues)
        k  = blues(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),2)  pos(tnp(k,2),2)], ...
             'Color', 'b', 'LineWidth', lw);
    end
end

% --- Black throats (fully liquid) ----------------------------------------
%  2D: if length(blacks)>0 → for z=1:length(blacks) → line(...,'color','k')
if ~isempty(blacks)
    for z = 1:length(blacks)
        k  = blacks(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),2)  pos(tnp(k,2),2)], ...
             'Color', 'k', 'LineWidth', lw);
    end
end

% --- White throats (fully vapour) ----------------------------------------
%  2D: if length(whites)>0 → for z=1:length(whites) → line(...,'color','w')
if ~isempty(whites)
    for z = 1:length(whites)
        k  = whites(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),2)  pos(tnp(k,2),2)], ...
             'Color', 'w', 'LineWidth', lw);
    end
end

% --- Title: same format as prof ------------------------------------------
%  2D: title(['Saturation = ',num2str(Ss(end),'%1.4f')],'FontSize',18)
%  Here Ss(s) is the current step value (Ss is being built up during loop).
title(['Saturation = ', num2str(Ss(s), '%1.4f'), '  (XY view)'], 'FontSize', 14);
hold off;
set(gcf, 'color', 'w');    % prof: set(gcf,'color','w')
axis off;
axis equal;

%% =========================================================================
%%  SUBPLOT 2  –  XZ PROJECTION  (top / depth view)
%%
%%  Axis: X = pos(:,1),  Z = pos(:,3)
%%  Same colour/width rules; box height = L4*(p-0.5) to match depth.
%%  No 2D analogue – added to expose the Z-direction drying front.
%% =========================================================================
subplot(1, 3, 2);

rectangle('Position', [-L4/2, -L4/2, L4*n, L4*(p-0.5)], ...
          'Curvature', [0 0], ...
          'FaceColor', [0.6 0.6 0.6], ...
          'EdgeColor', 'none');
hold on;

if ~isempty(blues)
    for z = 1:length(blues)
        k  = blues(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),3)  pos(tnp(k,2),3)], ...
             'Color', 'b', 'LineWidth', lw);
    end
end

if ~isempty(blacks)
    for z = 1:length(blacks)
        k  = blacks(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),3)  pos(tnp(k,2),3)], ...
             'Color', 'k', 'LineWidth', lw);
    end
end

if ~isempty(whites)
    for z = 1:length(whites)
        k  = whites(z);
        lw = lw_for(k);
        line([pos(tnp(k,1),1)  pos(tnp(k,2),1)], ...
             [pos(tnp(k,1),3)  pos(tnp(k,2),3)], ...
             'Color', 'w', 'LineWidth', lw);
    end
end

title(['Step = ', num2str(s), '  (XZ view)'], 'FontSize', 14);
hold off;
axis off;
axis equal;

%% =========================================================================
%%  SUBPLOT 3  –  LIVE DRYING CURVE  (exact replica of prof's subplot(1,2,2))
%%
%%  2D:  plot(Ss, Mev)
%%       title(['Saturation = ',num2str(Ss(end),'%1.4f')],'FontSize',18)
%%       axis([0 1 0 3e-10])
%%
%%  3D:  Identical.  Ss and Mev are built step-by-step in dryingLoop_3D;
%%       we plot only the non-zero portion (indices 1:s) so the curve
%%       grows in real-time exactly as in the 2D video.
%% =========================================================================
subplot(1, 3, 3);

% Plot only computed steps so far (same visual growth as 2D)
plot(Ss(1:s), Mev(1:s), 'b-', 'LineWidth', 1.2);

% Prof's exact title format
title(['Saturation = ', num2str(Ss(s), '%1.4f')], 'FontSize', 14);

% Prof's exact axis limits
axis([0 1 0 3e-10]);

xlabel('Saturation  S  [-]',        'FontSize', 11);
ylabel('Evap. rate  \dot{M}  [kg/s]','FontSize', 11);
grid on;
set(gca, 'XDir', 'reverse');   % drying goes right→left (same as 2D)