%% =========================================================================
%%  identifyCandidateThroats_3D.m
%%  Exact 3D replica of prof's candidate-throat identification file
%%
%%  PURPOSE
%%  -------
%%  Scans every throat in the network and:
%%    1. Computes Pore_Cond  – how many of the two endpoint pores are liquid
%%    2. Identifies invasion candidates  (liquid throat with ≥1 vapour pore)
%%    3. Sets Atactive, Pressure_meniscus for candidates
%%    4. Computes liquid conductance gw0 for all liquid-filled throats
%%    5. Returns MM  – index list of candidate throats (rcandidate ≠ 0)
%%
%%  REQUIRES (workspace)
%%  --------------------
%%  network.static : NT, tnp [NT×2], At [NT×1], Lt [NT×1], rt [Nt×1]
%%  ts(NT,1)       : throat state  (>0 = liquid-filled, 0 = vapour)
%%  ps(NP,1)       : pore  state  ( 1 = liquid,          0 = vapour)
%%  Patm           : atmospheric / reference pressure  [Pa]
%%  sigma          : surface tension                   [N/m]
%%  mu             : dynamic viscosity                 [Pa·s]
%%  RhoL           : liquid density                    [kg/m³]
%%
%%  PRODUCES
%%  --------
%%  Pore_Cond(NT,1)          : number of liquid endpoint pores (0, 1, or 2)
%%  Pore1(NT,1), Pore2(NT,1) : endpoint pore indices per throat
%%  rcandidate(NT,1)         : invasion-candidate radius (0 if not candidate)
%%  Atactive(NT,1)           : active cross-section area
%%  Pressure_meniscus(NT,1)  : capillary pressure at meniscus
%%  gw0(NT,1)                : liquid Hagen-Poiseuille conductance [kg/(s·Pa)]
%%  MM                       : indices of candidate throats
%% =========================================================================

%% --- Unpack ---------------------------------------------------------------
NT  = network.static.NT;
Nt  = network.static.Nt;
tnp = network.static.tnp;   % NT × 2
At  = network.static.At;    % NT × 1
Lt  = network.static.Lt;    % NT × 1
rt  = network.static.rt;    % Nt × 1  (interior throats only)

%% --- Reset output arrays -------------------------------------------------
%  (mirrors the initialisations already done in solveVaporPressure_3D;
%   called here again so this file is safe to call stand-alone mid-loop)
Pore_Cond         = zeros(NT, 1);
Pore1             = zeros(NT, 1);
Pore2             = zeros(NT, 1);
rcandidate        = zeros(NT, 1);
Atactive          = zeros(NT, 1);
Pressure_meniscus = zeros(NT, 1);
gw0               = zeros(NT, 1);

%% --- Main loop -----------------------------------------------------------
%  2D:  for i = 1:NT  (identical; NT now counts X+Y+Z+BL throats)
for i = 1:NT

    % ------------------------------------------------------------------
    % 2D: Pore_Cond(i) = sum(ps(nonzeros(tnp(i,:))))
    % tnp(i,:) has exactly 2 entries (never zero) so nonzeros is
    % redundant but kept for exact parity with the 2D formula.
    % ------------------------------------------------------------------
    Pore_Cond(i) = sum(ps(nonzeros(tnp(i,:))));

    % ------------------------------------------------------------------
    % 2D: Pore1(i) = tnp(i,1);  Pore2(i) = tnp(i,2);
    % UNCHANGED – tnp still has exactly 2 columns in 3D.
    % ------------------------------------------------------------------
    Pore1(i) = tnp(i,1);
    Pore2(i) = tnp(i,2);

    % ------------------------------------------------------------------
    % Invasion candidate:  liquid throat  AND  not both pores liquid
    % 2D: if ts(i) > 0 & (Pore_Cond(i) ~= 2)
    % UNCHANGED in logic.
    % ------------------------------------------------------------------
    if ts(i) > 0 && (Pore_Cond(i) ~= 2)
        rcandidate(i)        = rt(i);
        Atactive(i)          = At(i);
        Pressure_meniscus(i) = Patm - (2*sigma ./ rt(i));
    end

    % ------------------------------------------------------------------
    % Liquid conductance  (Hagen-Poiseuille, same formula as 2D)
    % 2D: if ts(i) ~= 0
    %       gw0(i) = (pi*rt(i)^4 / (mu*8*Lt(i)*ts(i))) * RhoL
    %     end
    % UNCHANGED. ts(i) is the liquid-saturation level of throat i;
    % dividing by ts(i) accounts for partial filling (same as 2D).
    % ------------------------------------------------------------------
    if ts(i) ~= 0
        gw0(i) = ((pi * rt(i).^4) ./ (mu * 8 * Lt(i) * ts(i))) * RhoL;
    end

end

%% --- Candidate index list ------------------------------------------------
% 2D: MM = find(rcandidate ~= 0)   UNCHANGED.
MM = find(rcandidate ~= 0);

fprintf('[identifyCandidateThroats_3D]  Candidates: %d / %d throats\n', ...
        length(MM), NT);