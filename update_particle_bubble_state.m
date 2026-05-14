function Particle = update_particle_bubble_state(Fluid, Particle, grid)
% =========================================================================
% update_particle_bubble_state
%
% Particle-bubble virtual contact state machine for unresolved liquid film.
%
% State definition:
%   0 = free
%   1 = collision
%   2 = induction
%   3 = attachment
%
% Main idea:
%   The real VOF field Fluid.alpha is not modified here.
%   This function only detects whether a particle is close enough to the
%   bubble interface and updates the particle-bubble interaction state.
%
% Distance definition used in this version:
%
%   h_alpha:
%       Minimum distance from real particle surface to current VOF interface
%       band, using the same idea as the visualization label:
%
%           alpha_low < alpha < alpha_high
%
%       By default:
%
%           alpha_low  = 0.05
%           alpha_high = 0.95
%
%   delta_v:
%       grid.virtual_film_thickness
%
%   h_p:
%       h_p = h_alpha - delta_v
%
% Important:
%   All state transitions in this function are based on h_p computed from
%   h_alpha, not from ray-searched h_geo.
%
%   Therefore the state criterion is consistent with the visualization label:
%
%       FREE:
%           h_p > h_c
%
%       COLLISION:
%           0 < h_p <= h_c
%
%       INDUCTION:
%           h_p <= 0
%
%       ATTACHMENT RELEASE:
%           h_p > grid.virtual_release_thickness
%
% Note:
%   This version intentionally keeps the original timer-reset behavior.
%   If a particle enters induction but in the next step h_p no longer
%   satisfies induction, the state is recomputed and the induction timer can
%   be reset to zero.
%
% Optional grid fields:
%   grid.enable_collision_state        default false
%   grid.enable_induction_state        default true
%   grid.virtual_film_thickness        default 1.5*h
%   grid.collision_hc_factor           default 0.2
%   grid.virtual_release_thickness     default 2.5*h
%   grid.pb_alpha_low                  default 0.05
%   grid.pb_alpha_high                 default 0.95
%   grid.induction_K                   default 1.0
%   grid.induction_t_min               default 2*dt
%   grid.induction_t_max               default 50*dt
%   grid.local_Rb_min                  default 2*h
%   grid.local_Rb_max_factor           default 1.0
%   grid.debug_pb_state                default false
%
% Particle fields initialized / updated:
%   Particle.pb_state
%   Particle.induction_timer
%   Particle.ta_induction
%   Particle.Rb_induction
%   Particle.virtual_contact_active
%   Particle.pb_contact_phi
%   Particle.pb_contact_x
%   Particle.pb_contact_y
%   Particle.pb_hp_geo       % kept for compatibility, here set equal to h_alpha
%   Particle.pb_hp_alpha     % actual alpha-band distance used by this function
%   Particle.pb_hp_model     % h_p = h_alpha - delta_v
%
% =========================================================================

%% ---------------- basic checks ----------------
if isempty(Particle) || ~isfield(Particle, 'r') || isempty(Particle.r)
    return;
end

if ~isfield(Fluid, 'alpha')
    error('update_particle_bubble_state: Fluid.alpha is required.');
end

alpha = Fluid.alpha;
alpha = apply_boundary_conditions('vof', grid, alpha);

Np = numel(Particle.r);

%% ---------------- initialize Particle fields ----------------
Particle = initialize_pb_state_fields(Particle, Np);

%% ---------------- grid / model parameters ----------------
h  = grid.h;
dt = grid.dt;

enable_collision_state = false;
if isfield(grid, 'enable_collision_state')
    enable_collision_state = logical(grid.enable_collision_state);
end

enable_induction_state = true;
if isfield(grid, 'enable_induction_state')
    enable_induction_state = logical(grid.enable_induction_state);
end

delta_v = 1.5 * h;
if isfield(grid, 'virtual_film_thickness')
    delta_v = grid.virtual_film_thickness;
end

hc_factor = 0.2;
if isfield(grid, 'collision_hc_factor')
    hc_factor = grid.collision_hc_factor;
end

release_thickness = 2.5 * h;
if isfield(grid, 'virtual_release_thickness')
    release_thickness = grid.virtual_release_thickness;
end

t_min = 2 * dt;
if isfield(grid, 'induction_t_min')
    t_min = grid.induction_t_min;
end

t_max = 50 * dt;
if isfield(grid, 'induction_t_max')
    t_max = grid.induction_t_max;
end

debug_pb_state = false;
if isfield(grid, 'debug_pb_state')
    debug_pb_state = logical(grid.debug_pb_state);
end

%% ---------------- loop over particles ----------------
for p = 1:Np

    old_state = Particle.pb_state(p);

    % ------------------------------------------------------------
    % 1. Find particle-bubble distance using alpha-band h_alpha
    % ------------------------------------------------------------
    gap_info = find_particle_bubble_gap_alpha(alpha, Particle, p, grid);

    if gap_info.found
        h_alpha = gap_info.hp_alpha;

        % State-machine distance:
        % This is the same definition as visualization:
        %
        %   h_p = gap_alpha - delta_v
        %
        h_p = h_alpha - delta_v;
    else
        h_alpha = inf;
        h_p     = inf;
    end

    % For compatibility with older post-processing:
    %   pb_hp_geo is no longer ray-searched h_geo here.
    %   It is set equal to h_alpha.
    Particle.pb_hp_geo(p)   = h_alpha;
    Particle.pb_hp_alpha(p) = h_alpha;
    Particle.pb_hp_model(p) = h_p;

    Particle.pb_contact_phi(p) = gap_info.phi;
    Particle.pb_contact_x(p)   = gap_info.x_if;
    Particle.pb_contact_y(p)   = gap_info.y_if;

    if debug_pb_state
        fprintf(['[PB state check] p=%d, old_state=%d, found=%d, ', ...
                 'h_alpha=%.3e m (%.2fh), h_p=%.3e m (%.2fh), ', ...
                 'delta_v=%.2fh\n'], ...
                 p, old_state, gap_info.found, ...
                 h_alpha, h_alpha/h, ...
                 h_p, h_p/h, ...
                 delta_v/h);
    end

    % ------------------------------------------------------------
    % 2. Attachment state has hysteresis / release condition
    % ------------------------------------------------------------
    if old_state == 3
        if ~gap_info.found || h_p > release_thickness
            Particle.pb_state(p) = 0;
            Particle.virtual_contact_active(p) = false;
            Particle.induction_timer(p) = 0.0;

            if debug_pb_state
                fprintf(['[PB release] p=%d, old_state=3 -> FREE, ', ...
                         'found=%d, h_p=%.3e m (%.2fh), ', ...
                         'release=%.3e m (%.2fh)\n'], ...
                         p, gap_info.found, ...
                         h_p, h_p/h, ...
                         release_thickness, release_thickness/h);
            end
        else
            Particle.pb_state(p) = 3;
            Particle.virtual_contact_active(p) = true;

            if debug_pb_state
                fprintf(['[PB keep attachment] p=%d, h_p=%.3e m (%.2fh), ', ...
                         'release=%.3e m (%.2fh)\n'], ...
                         p, h_p, h_p/h, ...
                         release_thickness, release_thickness/h);
            end
        end

        continue;
    end

    % ------------------------------------------------------------
    % 3. If no interface band found nearby, reset to free
    % ------------------------------------------------------------
    if ~gap_info.found
        Particle.pb_state(p) = 0;
        Particle.virtual_contact_active(p) = false;
        Particle.induction_timer(p) = 0.0;

        if debug_pb_state
            fprintf('[PB reset] p=%d, no alpha-band interface found -> FREE\n', p);
        end

        continue;
    end

    % ------------------------------------------------------------
    % 4. Determine candidate state from h_p based on h_alpha
    %
    %   free:      h_p > h_c
    %   collision: 0 < h_p <= h_c
    %   induction: h_p <= 0
    % ------------------------------------------------------------
    r_p = Particle.r(p);
    h_c = hc_factor * r_p;

    if h_p > h_c

        candidate_state = 0;   % free

    elseif h_p > 0 && h_p <= h_c

        if enable_collision_state
            candidate_state = 1;   % collision
        else
            candidate_state = 0;   % collision state disabled
        end

    else
        % h_p <= 0
        if enable_induction_state
            candidate_state = 2;   % induction
        else
            candidate_state = 3;   % directly attach for debugging
        end
    end

    if debug_pb_state
        fprintf(['[PB candidate] p=%d, h_p=%.3e m (%.2fh), ', ...
                 'h_c=%.3e m (%.2fh), candidate_state=%d, ', ...
                 'enable_collision=%d, enable_induction=%d\n'], ...
                 p, h_p, h_p/h, ...
                 h_c, h_c/h, ...
                 candidate_state, ...
                 enable_collision_state, enable_induction_state);
    end

    % ------------------------------------------------------------
    % 5. State transition logic
    %
    % Important:
    %   This keeps the original reset behavior.
    %   If old_state == 2 but current h_p no longer gives candidate_state=2,
    %   the timer will be reset in case 0 or case 1.
    % ------------------------------------------------------------
    switch candidate_state

        case 0
            % free
            Particle.pb_state(p) = 0;
            Particle.virtual_contact_active(p) = false;
            Particle.induction_timer(p) = 0.0;

            if debug_pb_state
                fprintf('[PB transition] p=%d -> FREE, timer reset\n', p);
            end

        case 1
            % collision
            Particle.pb_state(p) = 1;
            Particle.virtual_contact_active(p) = false;
            Particle.induction_timer(p) = 0.0;

            if debug_pb_state
                fprintf('[PB transition] p=%d -> COLLISION, timer reset\n', p);
            end

        case 2
            % induction
            if old_state ~= 2
                % Just entered induction. Freeze local bubble radius and ta.
                Particle.pb_state(p) = 2;
                Particle.virtual_contact_active(p) = false;
                Particle.induction_timer(p) = 0.0;

                Rb_loc = estimate_local_bubble_radius(alpha, gap_info, grid);
                Particle.Rb_induction(p) = Rb_loc;

                ta_raw = compute_induction_time(Rb_loc, r_p, grid);
                ta = min(max(ta_raw, t_min), t_max);
                Particle.ta_induction(p) = ta;

                if debug_pb_state
                    fprintf(['[PB enter induction] p=%d, h_p=%.3e m (%.2fh), ', ...
                             'Rb=%.3e m, rp=%.3e m, ', ...
                             'ta_raw=%.3e s, ta_clipped=%.3e s, ', ...
                             't_min=%.3e s, t_max=%.3e s, dt=%.3e s\n'], ...
                             p, h_p, h_p/h, ...
                             Rb_loc, r_p, ...
                             ta_raw, ta, ...
                             t_min, t_max, dt);
                end
            else
                Particle.pb_state(p) = 2;
                Particle.virtual_contact_active(p) = false;

                if debug_pb_state
                    fprintf(['[PB stay induction] p=%d, h_p=%.3e m (%.2fh), ', ...
                             'timer=%.3e s, ta=%.3e s\n'], ...
                             p, h_p, h_p/h, ...
                             Particle.induction_timer(p), ...
                             Particle.ta_induction(p));
                end
            end

            % Count induction time.
            Particle.induction_timer(p) = Particle.induction_timer(p) + dt;

            if debug_pb_state
                fprintf(['[PB induction timer] p=%d, h_p=%.3e m (%.2fh), ', ...
                         'timer=%.3e s / ta=%.3e s, dt=%.3e s\n'], ...
                         p, h_p, h_p/h, ...
                         Particle.induction_timer(p), ...
                         Particle.ta_induction(p), ...
                         dt);
            end

            if Particle.induction_timer(p) >= Particle.ta_induction(p)
                Particle.pb_state(p) = 3;
                Particle.virtual_contact_active(p) = true;

                if debug_pb_state
                    fprintf(['[PB attachment] p=%d attached. ', ...
                             'timer=%.3e s, ta=%.3e s, h_p=%.3e m (%.2fh)\n'], ...
                             p, ...
                             Particle.induction_timer(p), ...
                             Particle.ta_induction(p), ...
                             h_p, h_p/h);
                end
            end

        case 3
            % direct attachment, mainly used when induction is disabled
            Particle.pb_state(p) = 3;
            Particle.virtual_contact_active(p) = true;

            if Particle.Rb_induction(p) <= 0 || ~isfinite(Particle.Rb_induction(p))
                Particle.Rb_induction(p) = estimate_local_bubble_radius(alpha, gap_info, grid);
            end

            if Particle.ta_induction(p) <= 0 || ~isfinite(Particle.ta_induction(p))
                Particle.ta_induction(p) = compute_induction_time( ...
                    Particle.Rb_induction(p), r_p, grid);
            end

            if debug_pb_state
                fprintf(['[PB direct attachment] p=%d, h_p=%.3e m (%.2fh), ', ...
                         'Rb=%.3e m, ta=%.3e s\n'], ...
                         p, h_p, h_p/h, ...
                         Particle.Rb_induction(p), ...
                         Particle.ta_induction(p));
            end
    end
end

end

%% =========================================================================
% Subfunction: find_particle_bubble_gap_alpha
%% =========================================================================
function info = find_particle_bubble_gap_alpha(alpha, Particle, p, grid)
% =========================================================================
% find_particle_bubble_gap_alpha
%
% Find the closest VOF interface-band cell to the particle surface.
%
% This is designed to match the visualization gap definition:
%
%   interface band:
%       alpha_low < alpha < alpha_high
%
%   h_alpha:
%       min distance from particle surface to any interface-band cell center
%
%       h_alpha = min( distance(cell_center, particle_center) - r_p )
%
% This does NOT search alpha=0.5 along rays. It uses the same alpha-band
% standard as the label visualization.
%
% Output:
%   info.found
%   info.hp_alpha
%   info.hp_geo       % compatibility alias, set equal to hp_alpha
%   info.phi
%   info.normal
%   info.x_if, info.y_if
%   info.x_particle, info.y_particle
%   info.x_virtual, info.y_virtual
% =========================================================================

h        = grid.h;
start    = grid.start;
endy     = grid.endy;
endx     = grid.endx;
ghostnum = grid.ghostnum;

x_c = Particle.x_c(p);
y_c = Particle.y_c(p);
r_p = Particle.r(p);

delta_v = 1.5 * h;
if isfield(grid, 'virtual_film_thickness')
    delta_v = grid.virtual_film_thickness;
end

alpha_low = 0.05;
if isfield(grid, 'pb_alpha_low')
    alpha_low = grid.pb_alpha_low;
end

alpha_high = 0.95;
if isfield(grid, 'pb_alpha_high')
    alpha_high = grid.pb_alpha_high;
end

% Optional local search radius to avoid scanning whole domain too heavily.
% Keep large enough to cover near-field model range.
search_max_factor_h = 8.0;
if isfield(grid, 'pb_search_max_factor_h')
    search_max_factor_h = grid.pb_search_max_factor_h;
end

hc_factor = 0.2;
if isfield(grid, 'collision_hc_factor')
    hc_factor = grid.collision_hc_factor;
end

search_radius = r_p + max(search_max_factor_h*h, delta_v + hc_factor*r_p + 4*h);

% Initialize output.
info.found      = false;
info.hp_alpha   = inf;
info.hp_geo     = inf;
info.phi        = NaN;
info.normal     = [NaN; NaN];
info.x_if       = NaN;
info.y_if       = NaN;
info.x_particle = NaN;
info.y_particle = NaN;
info.x_virtual  = NaN;
info.y_virtual  = NaN;

best_gap = inf;

% Convert search box to index range.
ic = round(x_c/h + ghostnum + 0.5);
jc = round(y_c/h + ghostnum + 0.5);

Rbox = ceil(search_radius / h) + 2;

i1 = max(start, ic - Rbox);
i2 = min(endx,  ic + Rbox);
j1 = max(start, jc - Rbox);
j2 = min(endy,  jc + Rbox);

for j = j1:j2
    y = (j - ghostnum - 0.5) * h;

    for i = i1:i2
        a = alpha(j,i);

        if a <= alpha_low || a >= alpha_high
            continue;
        end

        x = (i - ghostnum - 0.5) * h;

        dx = x - x_c;
        dy = y - y_c;

        rr = hypot(dx, dy);

        % Ignore invalid exactly-centered case.
        if rr < 1e-30
            continue;
        end

        gap = rr - r_p;

        if gap < best_gap
            best_gap = gap;

            nx = dx / rr;
            ny = dy / rr;

            info.found    = true;
            info.hp_alpha = gap;
            info.hp_geo   = gap;

            info.phi    = atan2(ny, nx);
            info.normal = [nx; ny];

            info.x_if = x;
            info.y_if = y;

            info.x_particle = x_c + r_p * nx;
            info.y_particle = y_c + r_p * ny;

            info.x_virtual = x_c + (r_p + delta_v) * nx;
            info.y_virtual = y_c + (r_p + delta_v) * ny;
        end
    end
end

end

%% =========================================================================
% Subfunction: initialize fields
%% =========================================================================
function Particle = initialize_pb_state_fields(Particle, Np)

if ~isfield(Particle, 'pb_state') || numel(Particle.pb_state) ~= Np
    Particle.pb_state = zeros(Np,1);
else
    Particle.pb_state = Particle.pb_state(:);
end

if ~isfield(Particle, 'induction_timer') || numel(Particle.induction_timer) ~= Np
    Particle.induction_timer = zeros(Np,1);
else
    Particle.induction_timer = Particle.induction_timer(:);
end

if ~isfield(Particle, 'ta_induction') || numel(Particle.ta_induction) ~= Np
    Particle.ta_induction = zeros(Np,1);
else
    Particle.ta_induction = Particle.ta_induction(:);
end

if ~isfield(Particle, 'Rb_induction') || numel(Particle.Rb_induction) ~= Np
    Particle.Rb_induction = zeros(Np,1);
else
    Particle.Rb_induction = Particle.Rb_induction(:);
end

if ~isfield(Particle, 'virtual_contact_active') || numel(Particle.virtual_contact_active) ~= Np
    Particle.virtual_contact_active = false(Np,1);
else
    Particle.virtual_contact_active = logical(Particle.virtual_contact_active(:));
end

if ~isfield(Particle, 'pb_contact_phi') || numel(Particle.pb_contact_phi) ~= Np
    Particle.pb_contact_phi = nan(Np,1);
else
    Particle.pb_contact_phi = Particle.pb_contact_phi(:);
end

if ~isfield(Particle, 'pb_contact_x') || numel(Particle.pb_contact_x) ~= Np
    Particle.pb_contact_x = nan(Np,1);
else
    Particle.pb_contact_x = Particle.pb_contact_x(:);
end

if ~isfield(Particle, 'pb_contact_y') || numel(Particle.pb_contact_y) ~= Np
    Particle.pb_contact_y = nan(Np,1);
else
    Particle.pb_contact_y = Particle.pb_contact_y(:);
end

if ~isfield(Particle, 'pb_hp_geo') || numel(Particle.pb_hp_geo) ~= Np
    Particle.pb_hp_geo = inf(Np,1);
else
    Particle.pb_hp_geo = Particle.pb_hp_geo(:);
end

if ~isfield(Particle, 'pb_hp_alpha') || numel(Particle.pb_hp_alpha) ~= Np
    Particle.pb_hp_alpha = inf(Np,1);
else
    Particle.pb_hp_alpha = Particle.pb_hp_alpha(:);
end

if ~isfield(Particle, 'pb_hp_model') || numel(Particle.pb_hp_model) ~= Np
    Particle.pb_hp_model = inf(Np,1);
else
    Particle.pb_hp_model = Particle.pb_hp_model(:);
end

end

%% =========================================================================
% Subfunction: estimate local bubble radius
%% =========================================================================
function Rb_loc = estimate_local_bubble_radius(alpha, gap_info, grid)
% =========================================================================
% estimate_local_bubble_radius
%
% Estimate local bubble radius near the detected particle-bubble contact
% region using local curvature from the real alpha field.
%
% Rb_loc = 1 / |kappa_loc|
%
% =========================================================================

h        = grid.h;
start    = grid.start;
endy     = grid.endy;
endx     = grid.endx;
Ny       = grid.Ny;
Nx       = grid.Nx;
ghostnum = grid.ghostnum;

x0 = gap_info.x_if;
y0 = gap_info.y_if;

if ~isfinite(x0) || ~isfinite(y0)
    Rb_loc = estimate_equivalent_bubble_radius(alpha, grid);
    return;
end

% Smoothing alpha slightly for curvature estimation.
alpha_s = alpha;
for j = start:endy
    for i = start:endx
        alpha_s(j,i) = ( ...
            1*alpha(j-1,i-1) + 2*alpha(j-1,i) + 1*alpha(j-1,i+1) + ...
            2*alpha(j  ,i-1) + 4*alpha(j  ,i) + 2*alpha(j  ,i+1) + ...
            1*alpha(j+1,i-1) + 2*alpha(j+1,i) + 1*alpha(j+1,i+1) ) / 16;
    end
end
alpha_s = apply_boundary_conditions('vof', grid, alpha_s);

nx = zeros(Ny, Nx);
ny = zeros(Ny, Nx);

for j = start:endy
    for i = start:endx
        ax = (alpha_s(j,i+1) - alpha_s(j,i-1)) / (2*h);
        ay = (alpha_s(j+1,i) - alpha_s(j-1,i)) / (2*h);

        gmag = sqrt(ax^2 + ay^2) + 1e-20;

        nx(j,i) = ax / gmag;
        ny(j,i) = ay / gmag;
    end
end

kappa = zeros(Ny, Nx);
for j = start+1:endy-1
    for i = start+1:endx-1
        kappa(j,i) = - ( ...
            (nx(j,i+1) - nx(j,i-1)) / (2*h) + ...
            (ny(j+1,i) - ny(j-1,i)) / (2*h) );
    end
end

% Weighted local averaging around the contact point.
R_window = 4.0*h;
if isfield(grid, 'local_Rb_window')
    R_window = grid.local_Rb_window;
end

sigma_w = 2.0*h;

ic = round(x0/h + ghostnum + 0.5);
jc = round(y0/h + ghostnum + 0.5);

Rbox = ceil(R_window/h) + 2;

i1 = max(start+1, ic - Rbox);
i2 = min(endx-1,  ic + Rbox);
j1 = max(start+1, jc - Rbox);
j2 = min(endy-1,  jc + Rbox);

num = 0.0;
den = 0.0;

for j = j1:j2
    for i = i1:i2

        a = alpha_s(j,i);

        % Interface cells only.
        if a <= 0.05 || a >= 0.95
            continue;
        end

        xc = (i - ghostnum - 0.5)*h;
        yc = (j - ghostnum - 0.5)*h;

        dist2 = (xc - x0)^2 + (yc - y0)^2;

        if dist2 > R_window^2
            continue;
        end

        w_if = a * (1 - a);
        w_r  = exp(-dist2 / (2*sigma_w^2));

        w = w_if * w_r;

        if isfinite(kappa(j,i))
            num = num + w * kappa(j,i);
            den = den + w;
        end
    end
end

if den <= 1e-30
    Rb_loc = estimate_equivalent_bubble_radius(alpha, grid);
else
    kappa_loc = num / den;

    if abs(kappa_loc) < 1e-12 || ~isfinite(kappa_loc)
        Rb_loc = estimate_equivalent_bubble_radius(alpha, grid);
    else
        Rb_loc = 1.0 / abs(kappa_loc);
    end
end

% Radius clipping.
Rb_min = 2.0*h;
if isfield(grid, 'local_Rb_min')
    Rb_min = grid.local_Rb_min;
end

Rb_equiv = estimate_equivalent_bubble_radius(alpha, grid);

Rb_max_factor = 1.0;
if isfield(grid, 'local_Rb_max_factor')
    Rb_max_factor = grid.local_Rb_max_factor;
end

Rb_max = max(Rb_min, Rb_max_factor * Rb_equiv);

Rb_loc = min(max(Rb_loc, Rb_min), Rb_max);

end

%% =========================================================================
% Subfunction: compute induction time
%% =========================================================================
function ta = compute_induction_time(Rb_loc, rp, grid)
% =========================================================================
% compute_induction_time
%
% Literature-style empirical induction time:
%
%   ta = K * db^2.38 * dp^1.59 / (db + dp)^1.59
%
% where:
%   db = 2*Rb_loc
%   dp = 2*rp
%
% K must be calibrated in the user's unit system.
% =========================================================================

K = 1.0;
if isfield(grid, 'induction_K')
    K = grid.induction_K;
end

db = 2.0 * Rb_loc;
dp = 2.0 * rp;

db = max(db, 1e-30);
dp = max(dp, 1e-30);

ta = K * (db^2.38 * dp^1.59) / ((db + dp)^1.59 + 1e-30);

if ~isfinite(ta)
    ta = 0.0;
end

end

%% =========================================================================
% Subfunction: equivalent bubble radius
%% =========================================================================
function Rb_equiv = estimate_equivalent_bubble_radius(alpha, grid)
% Estimate an equivalent bubble radius from gas area.
% Assumption:
%   alpha = 1 liquid, alpha = 0 gas.
%
% For multiple bubbles this gives a global effective radius. This is only a
% fallback. For multiple bubbles, a connected-component bubble radius should
% be implemented later.

h     = grid.h;
start = grid.start;
endy  = grid.endy;
endx  = grid.endx;

alpha_active = alpha(start:endy, start:endx);

gas_fraction = 1.0 - alpha_active;
gas_fraction = min(max(gas_fraction, 0.0), 1.0);

A_gas = sum(gas_fraction(:)) * h^2;

if A_gas <= 1e-30
    Rb_equiv = 4*h;
else
    Rb_equiv = sqrt(A_gas / pi);
end

Rb_min = 2*h;
if isfield(grid, 'local_Rb_min')
    Rb_min = grid.local_Rb_min;
end

Rb_equiv = max(Rb_equiv, Rb_min);

end
