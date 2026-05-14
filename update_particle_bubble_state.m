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
% Definitions:
%   h_geo   = geometrical gap between real particle surface and alpha=0.5
%             bubble interface, searched along particle outward normals.
%
%   delta_v = grid.virtual_film_thickness
%
%   h_p_model = h_geo - delta_v
%
% Interpretation:
%   h_p_model < 0 means the real geometrical gap is smaller than the
%   unresolved virtual film thickness. This is used as the numerical
%   counterpart of the h_p < 0 criterion in the induction-stage model.
%
% Required / suggested grid fields:
%   grid.h
%   grid.dt
%   grid.Nx, grid.Ny
%   grid.ghostnum
%   grid.start, grid.endx, grid.endy
%
% Optional grid fields:
%   grid.enable_collision_state        default true
%   grid.enable_induction_state        default true
%   grid.virtual_film_thickness        default 1.5*h
%   grid.collision_hc_factor           default 0.2
%   grid.virtual_release_thickness     default 2.5*h
%   grid.pb_N_lag                      default 360
%   grid.pb_search_ds_factor           default 0.25
%   grid.pb_search_max_factor_h        default 8
%   grid.pb_alpha_interface            default 0.5
%   grid.induction_K                   default 1.0
%   grid.induction_t_min               default 2*dt
%   grid.induction_t_max               default 50*dt
%   grid.local_Rb_min                  default 2*h
%   grid.local_Rb_max_factor           default 1.0
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
%   Particle.pb_hp_geo
%   Particle.pb_hp_model
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

%% ---------------- loop over particles ----------------
for p = 1:Np

    old_state = Particle.pb_state(p);

    % ------------------------------------------------------------
    % 1. Find geometrical particle-bubble gap using real alpha field
    % ------------------------------------------------------------
    gap_info = find_particle_bubble_gap(alpha, Particle, p, grid);

    if gap_info.found
        h_geo   = gap_info.hp_geo;
        h_model = h_geo - delta_v;
    else
        h_geo   = inf;
        h_model = inf;
    end

    Particle.pb_hp_geo(p)   = h_geo;
    Particle.pb_hp_model(p) = h_model;

    Particle.pb_contact_phi(p) = gap_info.phi;
    Particle.pb_contact_x(p)   = gap_info.x_if;
    Particle.pb_contact_y(p)   = gap_info.y_if;

    % ------------------------------------------------------------
    % 2. Attachment state has hysteresis / release condition
    % ------------------------------------------------------------
    if old_state == 3
        if ~gap_info.found || h_model > release_thickness
            Particle.pb_state(p) = 0;
            Particle.virtual_contact_active(p) = false;
            Particle.induction_timer(p) = 0.0;
        else
            Particle.pb_state(p) = 3;
            Particle.virtual_contact_active(p) = true;
        end

        continue;
    end

    % ------------------------------------------------------------
    % 3. If no interface found nearby, reset to free
    % ------------------------------------------------------------
    if ~gap_info.found
        Particle.pb_state(p) = 0;
        Particle.virtual_contact_active(p) = false;
        Particle.induction_timer(p) = 0.0;
        continue;
    end

    % ------------------------------------------------------------
    % 4. Determine candidate state from h_p_model
    % ------------------------------------------------------------
    r_p = Particle.r(p);
    h_c = hc_factor * r_p;

    if h_model > h_c

        candidate_state = 0;   % free

    elseif h_model > 0 && h_model <= h_c

        if enable_collision_state
            candidate_state = 1;   % collision
        else
            candidate_state = 0;   % collision state disabled
        end

    else
        % h_model <= 0
        if enable_induction_state
            candidate_state = 2;   % induction
        else
            candidate_state = 3;   % directly attach for debugging
        end
    end

    % ------------------------------------------------------------
    % 5. State transition logic
    % ------------------------------------------------------------
    switch candidate_state

        case 0
            % free
            Particle.pb_state(p) = 0;
            Particle.virtual_contact_active(p) = false;
            Particle.induction_timer(p) = 0.0;

        case 1
            % collision
            Particle.pb_state(p) = 1;
            Particle.virtual_contact_active(p) = false;

            % Keep timer at zero before induction.
            Particle.induction_timer(p) = 0.0;

        case 2
            % induction
            if old_state ~= 2
                % Just entered induction. Freeze local bubble radius and ta.
                Particle.pb_state(p) = 2;
                Particle.virtual_contact_active(p) = false;
                Particle.induction_timer(p) = 0.0;

                Rb_loc = estimate_local_bubble_radius(alpha, gap_info, grid);
                Particle.Rb_induction(p) = Rb_loc;

                ta = compute_induction_time(Rb_loc, r_p, grid);
                ta = min(max(ta, t_min), t_max);
                Particle.ta_induction(p) = ta;
            else
                Particle.pb_state(p) = 2;
                Particle.virtual_contact_active(p) = false;
            end

            % Count induction time.
            Particle.induction_timer(p) = Particle.induction_timer(p) + dt;

            if Particle.induction_timer(p) >= Particle.ta_induction(p)
                Particle.pb_state(p) = 3;
                Particle.virtual_contact_active(p) = true;
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
    end
end

end

%% =========================================================================
% Subfunction: find_particle_bubble_gap
%% =========================================================================
function info = find_particle_bubble_gap(alpha, Particle, p, grid)
% =========================================================================
% find_particle_bubble_gap
%
% Search alpha=0.5 interface along outward rays from particle surface.
%
% For each angular direction:
%   x(s) = x_c + (r_p + s) * n
%   y(s) = y_c + (r_p + s) * n
%
% where s >= 0 is the distance from the real particle surface.
%
% The nearest alpha=0.5 crossing is returned as h_geo.
%
% Output info:
%   info.found
%   info.hp_geo
%   info.phi
%   info.normal
%   info.x_if, info.y_if
%   info.x_particle, info.y_particle
%   info.x_virtual, info.y_virtual
% =========================================================================

h        = grid.h;
Ny       = grid.Ny;
Nx       = grid.Nx;
ghostnum = grid.ghostnum;

x_c = Particle.x_c(p);
y_c = Particle.y_c(p);
r_p = Particle.r(p);

N_lag = 360;
if isfield(grid, 'pb_N_lag')
    N_lag = grid.pb_N_lag;
end

alpha_if = 0.5;
if isfield(grid, 'pb_alpha_interface')
    alpha_if = grid.pb_alpha_interface;
end

delta_v = 1.5 * h;
if isfield(grid, 'virtual_film_thickness')
    delta_v = grid.virtual_film_thickness;
end

ds_factor = 0.25;
if isfield(grid, 'pb_search_ds_factor')
    ds_factor = grid.pb_search_ds_factor;
end
ds = ds_factor * h;

search_max_factor_h = 8.0;
if isfield(grid, 'pb_search_max_factor_h')
    search_max_factor_h = grid.pb_search_max_factor_h;
end

hc_factor = 0.2;
if isfield(grid, 'collision_hc_factor')
    hc_factor = grid.collision_hc_factor;
end

% Search distance should cover virtual film + collision range + buffer.
s_max = max(search_max_factor_h * h, delta_v + hc_factor * r_p + 4*h);

N_search = max(2, ceil(s_max / ds));

% Initialize output.
info.found      = false;
info.hp_geo     = inf;
info.phi        = NaN;
info.normal     = [NaN; NaN];
info.x_if       = NaN;
info.y_if       = NaN;
info.x_particle = NaN;
info.y_particle = NaN;
info.x_virtual  = NaN;
info.y_virtual  = NaN;

best_s = inf;

phi_list = linspace(0, 2*pi, N_lag+1);
phi_list = phi_list(1:end-1);

for k = 1:N_lag

    phi = phi_list(k);
    nvec = [cos(phi); sin(phi)];

    % Values along the outward ray.
    s_prev = 0.0;
    x_prev = x_c + (r_p + s_prev) * nvec(1);
    y_prev = y_c + (r_p + s_prev) * nvec(2);

    a_prev = interp_cell_scalar_bilinear_local( ...
        alpha, x_prev, y_prev, h, ghostnum, Ny, Nx);

    % March outward.
    for m = 1:N_search

        s_curr = m * ds;

        x_curr = x_c + (r_p + s_curr) * nvec(1);
        y_curr = y_c + (r_p + s_curr) * nvec(2);

        % If outside active computational domain, stop this ray.
        if ~point_inside_active_domain(x_curr, y_curr, grid)
            break;
        end

        a_curr = interp_cell_scalar_bilinear_local( ...
            alpha, x_curr, y_curr, h, ghostnum, Ny, Nx);

        f_prev = a_prev - alpha_if;
        f_curr = a_curr - alpha_if;

        % Detect crossing of alpha=0.5.
        crossed = false;

        if f_prev == 0
            crossed = true;
            lambda = 0.0;
        elseif f_prev * f_curr < 0
            crossed = true;
            lambda = abs(f_prev) / (abs(f_prev) + abs(f_curr) + 1e-30);
        elseif f_curr == 0
            crossed = true;
            lambda = 1.0;
        else
            lambda = NaN;
        end

        if crossed
            s_if = (1 - lambda) * s_prev + lambda * s_curr;

            if s_if < best_s
                best_s = s_if;

                x_if = x_c + (r_p + s_if) * nvec(1);
                y_if = y_c + (r_p + s_if) * nvec(2);

                info.found      = true;
                info.hp_geo     = s_if;
                info.phi        = phi;
                info.normal     = nvec;
                info.x_if       = x_if;
                info.y_if       = y_if;
                info.x_particle = x_c + r_p * nvec(1);
                info.y_particle = y_c + r_p * nvec(2);
                info.x_virtual  = x_c + (r_p + delta_v) * nvec(1);
                info.y_virtual  = y_c + (r_p + delta_v) * nvec(2);
            end

            % The first crossing along this ray is enough.
            break;
        end

        s_prev = s_curr;
        a_prev = a_curr;
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
% This is intentionally a simple and robust first version. The radius is
% frozen at the moment of entering induction in the main state function.
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

%% =========================================================================
% Subfunction: bilinear interpolation of cell-centered scalar
%% =========================================================================
function val = interp_cell_scalar_bilinear_local(phi, xq, yq, h, ghostnum, Ny, Nx)

i_idx = xq / h + 0.5 + ghostnum;
j_idx = yq / h + 0.5 + ghostnum;

i0 = floor(i_idx);
j0 = floor(j_idx);
i1 = i0 + 1;
j1 = j0 + 1;

i0 = max(1, min(i0, Nx));
i1 = max(1, min(i1, Nx));
j0 = max(1, min(j0, Ny));
j1 = max(1, min(j1, Ny));

dx = i_idx - i0;
dy = j_idx - j0;

dx = min(max(dx, 0.0), 1.0);
dy = min(max(dy, 0.0), 1.0);

val = (1-dx)*(1-dy)*phi(j0,i0) + ...
       dx   *(1-dy)*phi(j0,i1) + ...
      (1-dx)* dy   *phi(j1,i0) + ...
       dx   * dy   *phi(j1,i1);

val = min(max(val, 0.0), 1.0);

end

%% =========================================================================
% Subfunction: active-domain check
%% =========================================================================
function tf = point_inside_active_domain(x, y, grid)

h        = grid.h;
ghostnum = grid.ghostnum;
start    = grid.start;
endy     = grid.endy;
endx     = grid.endx;

x_min = (start - ghostnum - 0.5)*h;
x_max = (endx  - ghostnum - 0.5)*h;
y_min = (start - ghostnum - 0.5)*h;
y_max = (endy  - ghostnum - 0.5)*h;

tf = (x >= x_min && x <= x_max && y >= y_min && y <= y_max);

end