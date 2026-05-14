function alpha_ext = update_extension(alpha_in, Particle, grid)
% =========================================================================
% update_extension
%
% Virtual-contact Liu-style characteristic extension for unresolved
% particle-bubble liquid film.
%
% This version is modified from the original update_extension.m.
%
% Key differences from the original version:
%   1) This function now only works for particles with
%          Particle.virtual_contact_active(p) == true
%
%   2) Contact-line detection is performed on a virtual particle interface:
%          r_detect = r_real + grid.virtual_film_thickness
%
%      This corresponds to unresolved-film / induction-attachment modeling.
%
%   3) The extension is still written only inside the real solid particle:
%          rr < r_real
%
%      Therefore alpha_ext is only a capillary-geometry auxiliary field and
%      does not overwrite the real VOF field Fluid.alpha.
%
%   4) The internal extension depth is still controlled by:
%          grid.extension_distance
%
% Recommended usage in the main loop:
%   alpha_real_new = update_advection(Fluid, grid);
%   alpha_ext_new  = update_extension(alpha_real_new, Particle, grid);
%
% Then:
%   Fluid.alpha     = alpha_real_new;
%   Fluid.alpha_ext = alpha_ext_new;
%
% prediction.m should use Fluid.alpha_ext only for CSF if desired.
% Density / viscosity / mass conservation should still use Fluid.alpha.
% =========================================================================

Ny       = grid.Ny;
Nx       = grid.Nx;
h        = grid.h;
ghostnum = grid.ghostnum;
start    = grid.start;
endy     = grid.endy;
endx     = grid.endx;

alpha_ext = alpha_in;
alpha_ext = apply_boundary_conditions('vof', grid, alpha_ext);

if isempty(Particle) || ~isfield(Particle,'r') || isempty(Particle.r)
    return;
end

% -------------------------------------------------------------------------
% Global virtual-extension switch
% -------------------------------------------------------------------------
enable_virtual_extension = true;
if isfield(grid, 'enable_virtual_extension')
    enable_virtual_extension = logical(grid.enable_virtual_extension);
end

if ~enable_virtual_extension
    return;
end

% This function is now state-machine driven.
% If virtual_contact_active does not exist, do nothing.
if ~isfield(Particle, 'virtual_contact_active')
    return;
end

if ~any(Particle.virtual_contact_active)
    return;
end

eps_div = 1e-14;

% -------------------------------------------------------------------------
% alpha convention
% -------------------------------------------------------------------------
alpha_represents_liquid = true;
if isfield(grid,'alpha_represents_liquid')
    alpha_represents_liquid = logical(grid.alpha_represents_liquid);
end

% -------------------------------------------------------------------------
% interface detector
% -------------------------------------------------------------------------
alpha_min = 1e-3;
alpha_max = 1 - 1e-3;

if isfield(grid,'capillary_alpha_min')
    alpha_min = grid.capillary_alpha_min;
end

if isfield(grid,'capillary_alpha_max')
    alpha_max = grid.capillary_alpha_max;
end

% -------------------------------------------------------------------------
% extension parameters
% -------------------------------------------------------------------------
ext_dist = 2.0*h;
if isfield(grid,'extension_distance')
    ext_dist = grid.extension_distance;
end

band_tan_halfwidth = 3.0*h;
if isfield(grid,'extension_tangent_halfwidth')
    band_tan_halfwidth = grid.extension_tangent_halfwidth;
end

N_lag_default = 240;
if isfield(grid,'capillary_N_lag')
    N_lag_default = grid.capillary_N_lag;
end

probe_dist = 0.5*h;
if isfield(grid,'capillary_probe_dist')
    probe_dist = grid.capillary_probe_dist;
end

virtual_film_thickness = 1.5*h;
if isfield(grid, 'virtual_film_thickness')
    virtual_film_thickness = grid.virtual_film_thickness;
end

% If true, still require nearby interfacial alpha around the extension cell.
% For virtual contact this is usually too strict, because the real solid
% boundary and the virtual detection boundary are separated by delta_v.
require_interfacial_neighbor = false;
if isfield(grid, 'extension_require_interfacial_neighbor')
    require_interfacial_neighbor = logical(grid.extension_require_interfacial_neighbor);
end

% Optional fallback.
allow_bilinear_fallback = false;
if isfield(grid, 'extension_allow_bilinear_fallback')
    allow_bilinear_fallback = logical(grid.extension_allow_bilinear_fallback);
end

% Optional cleanup of deep solid residual alpha.
clean_deep_solid_alpha = false;
if isfield(grid, 'clean_deep_solid_alpha')
    clean_deep_solid_alpha = logical(grid.clean_deep_solid_alpha);
end

% Keep nearest contact-point contribution.
best_metric     = inf(Ny, Nx);
alpha_candidate = alpha_in;

Np = length(Particle.r);

for p = 1:Np

    % ---------------------------------------------------------------------
    % Only attachment-state particles do virtual extension.
    % ---------------------------------------------------------------------
    if ~Particle.virtual_contact_active(p)
        continue;
    end

    x_c = Particle.x_c(p);
    y_c = Particle.y_c(p);

    r_real   = Particle.r(p);
    r_detect = r_real + virtual_film_thickness;

    theta_eq = get_particle_contact_angle_local(Particle, p);

    % ---------------------------------------------------------------------
    % Lagrangian points on virtual circular immersed boundary
    % ---------------------------------------------------------------------
    N_lag = max(N_lag_default, ceil(2*pi*r_detect/(0.5*h)));

    ang = linspace(0, 2*pi, N_lag+1);
    ang = ang(1:end-1);

    x_lag = x_c + r_detect*cos(ang);
    y_lag = y_c + r_detect*sin(ang);

    alpha_probe = zeros(1, N_lag);

    for k = 1:N_lag
        xk = x_lag(k);
        yk = y_lag(k);

        ns = [xk - x_c; yk - y_c];
        ns = ns / (norm(ns) + eps_div);

        xq = xk + probe_dist*ns(1);
        yq = yk + probe_dist*ns(2);

        alpha_probe(k) = interp_cell_scalar_bilinear_local( ...
            alpha_in, xq, yq, h, ghostnum, Ny, Nx);
    end

    % ---------------------------------------------------------------------
    % Detect alpha = 0.5 crossings along the virtual boundary sampling.
    % ---------------------------------------------------------------------
    for k = 1:N_lag

        kp = k + 1;
        if kp > N_lag
            kp = 1;
        end

        s1 = alpha_probe(k)  - 0.5;
        s2 = alpha_probe(kp) - 0.5;

        near_contact = ...
            ((alpha_probe(k)  > alpha_min && alpha_probe(k)  < alpha_max) || ...
             (alpha_probe(kp) > alpha_min && alpha_probe(kp) < alpha_max));

        if ~near_contact
            continue;
        end

        if s1*s2 > 0 || abs(s1 - s2) < 1e-12
            continue;
        end

        lam = abs(s1) / (abs(s1) + abs(s2) + eps_div);

        % Contact point on virtual boundary.
        x_cp = (1-lam)*x_lag(k) + lam*x_lag(kp);
        y_cp = (1-lam)*y_lag(k) + lam*y_lag(kp);

        Xcp = [x_cp; y_cp];

        ns = [x_cp - x_c; y_cp - y_c];
        ns = ns / (norm(ns) + eps_div);

        % Tangent on the virtual particle surface.
        tw = [-ns(2); ns(1)];

        % -----------------------------------------------------------------
        % Effective geometric angle
        % -----------------------------------------------------------------
        if alpha_represents_liquid
            theta_geom = theta_eq;
        else
            theta_geom = pi - theta_eq;
        end

        % Two characteristic directions symmetric about the wall normal.
        m1 = rotate_vec_local(tw,  theta_geom);
        m2 = rotate_vec_local(tw, -theta_geom);

        m1 = m1 / (norm(m1) + eps_div);
        m2 = m2 / (norm(m2) + eps_div);

        % Local box around virtual contact point.
        ic = round(x_cp/h + ghostnum + 0.5);
        jc = round(y_cp/h + ghostnum + 0.5);

        Rbox = ceil((virtual_film_thickness + ext_dist + band_tan_halfwidth + 3*h)/h) + 4;

        i1 = max(start, ic - Rbox);
        i2 = min(endx,  ic + Rbox);
        j1 = max(start, jc - Rbox);
        j2 = min(endy,  jc + Rbox);

        for j = j1:j2
            for i = i1:i2

                xc = (i - ghostnum - 0.5)*h;
                yc = (j - ghostnum - 0.5)*h;
                P  = [xc; yc];

                vec_pc = P - [x_c; y_c];
                rr = norm(vec_pc);

                % ---------------------------------------------------------
                % Important:
                %   Only write extension into the real solid particle.
                %   Do not write into the virtual layer outside the particle.
                % ---------------------------------------------------------
                if rr >= r_real
                    continue;
                end

                % Only within real-solid extension band.
                d_ib_real = r_real - rr;
                if d_ib_real > ext_dist
                    continue;
                end

                rel = P - Xcp;

                s_tan = dot(rel, tw);

                % Into virtual solid positive.
                % Since Xcp is on the virtual surface, real solid cells will
                % have positive s_nrm approximately equal to delta_v + depth.
                s_nrm_virtual = -dot(rel, ns);

                if s_nrm_virtual < -0.5*h
                    continue;
                end

                if abs(s_tan) > band_tan_halfwidth
                    continue;
                end

                if require_interfacial_neighbor
                    if ~is_ghost_contact_band_cell(i, j, alpha_in, x_c, y_c, r_detect, ...
                            alpha_min, alpha_max, h, ghostnum, Ny, Nx)
                        continue;
                    end
                end

                % ---------------------------------------------------------
                % Characteristic interpolation based on the virtual surface.
                %
                % P is inside the virtual circle. D1/D2 are chosen outside
                % the virtual circle, so the interpolated alpha values are
                % taken from the real fluid side outside the unresolved film.
                % ---------------------------------------------------------
                Cc = [x_c; y_c];

                [ok1, CD1] = characteristic_quadric_value_circle( ...
                    alpha_in, P, m1, Cc, r_detect, h, ghostnum, Ny, Nx, start, endx, endy);

                [ok2, CD2] = characteristic_quadric_value_circle( ...
                    alpha_in, P, m2, Cc, r_detect, h, ghostnum, Ny, Nx, start, endx, endy);

                if (~ok1 || ~ok2) && allow_bilinear_fallback
                    [ok1b, D1b] = first_outside_point_simple(P, m1, Cc, r_detect);
                    [ok2b, D2b] = first_outside_point_simple(P, m2, Cc, r_detect);

                    if ~ok1 && ok1b
                        CD1 = interp_cell_scalar_bilinear_local(alpha_in, D1b(1), D1b(2), h, ghostnum, Ny, Nx);
                        ok1 = true;
                    end

                    if ~ok2 && ok2b
                        CD2 = interp_cell_scalar_bilinear_local(alpha_in, D2b(1), D2b(2), h, ghostnum, Ny, Nx);
                        ok2 = true;
                    end
                end

                if ~ok1 || ~ok2
                    continue;
                end

                % ---------------------------------------------------------
                % Liu's max/min rule
                % ---------------------------------------------------------
                if alpha_represents_liquid
                    if theta_eq <= pi/2
                        a_ext = max(CD1, CD2);
                    else
                        a_ext = min(CD1, CD2);
                    end
                else
                    if theta_eq <= pi/2
                        a_ext = min(CD1, CD2);
                    else
                        a_ext = max(CD1, CD2);
                    end
                end

                a_ext = min(max(a_ext, 0.0), 1.0);

                % Prefer closer to current contact point and to real solid boundary.
                metric = d_ib_real + 0.2*abs(s_tan);

                if metric < best_metric(j,i)
                    best_metric(j,i)     = metric;
                    alpha_candidate(j,i) = a_ext;
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % Optional: clean deep real-solid residual alpha for this particle.
    % ---------------------------------------------------------------------
    if clean_deep_solid_alpha
        for j = start:endy
            for i = start:endx

                xc = (i - ghostnum - 0.5)*h;
                yc = (j - ghostnum - 0.5)*h;

                rr = hypot(xc - x_c, yc - y_c);

                if rr < r_real
                    d_ib_real = r_real - rr;

                    if d_ib_real > ext_dist
                        if alpha_represents_liquid
                            alpha_candidate(j,i) = 0.0;
                        else
                            alpha_candidate(j,i) = 1.0;
                        end
                    end
                end
            end
        end
    end
end

alpha_ext = alpha_candidate;
alpha_ext(start:endy, start:endx) = min(max(alpha_ext(start:endy, start:endx), 0.0), 1.0);
alpha_ext = apply_boundary_conditions('vof', grid, alpha_ext);

end

%% =========================================================================
% helpers
%% =========================================================================

function theta_eq = get_particle_contact_angle_local(Particle, p)

if isfield(Particle, 'theta')
    theta_eq = Particle.theta(p);
elseif isfield(Particle, 'contact_angle')
    theta_eq = Particle.contact_angle(p);
elseif isfield(Particle, 'theta_eq')
    theta_eq = Particle.theta_eq(p);
else
    theta_eq = pi/2;
end

if theta_eq > 2*pi
    theta_eq = theta_eq * pi / 180.0;
end

end

function vout = rotate_vec_local(vin, ang)

c = cos(ang);
s = sin(ang);

R = [c, -s; ...
     s,  c];

vout = R * vin(:);

end

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

function flag = is_ghost_contact_band_cell(i, j, alpha, xc, yc, r, ...
                                           alpha_min, alpha_max, h, ghostnum, Ny, Nx)
% Rough analogue of Liu ghost contact-line region:
% point is inside solid/virtual-solid and itself or outside-neighbor is interfacial.

flag = false;

x0 = (i - ghostnum - 0.5)*h;
y0 = (j - ghostnum - 0.5)*h;
rr0 = hypot(x0 - xc, y0 - yc);

if rr0 >= r
    return;
end

if alpha(j,i) > alpha_min && alpha(j,i) < alpha_max
    flag = true;
    return;
end

nbr = [-1 0; ...
        1 0; ...
        0 -1; ...
        0 1];

for q = 1:size(nbr,1)

    ii = i + nbr(q,1);
    jj = j + nbr(q,2);

    if ii < 1 || ii > Nx || jj < 1 || jj > Ny
        continue;
    end

    x1 = (ii - ghostnum - 0.5)*h;
    y1 = (jj - ghostnum - 0.5)*h;
    rr1 = hypot(x1 - xc, y1 - yc);

    if rr1 >= r && alpha(jj,ii) > alpha_min && alpha(jj,ii) < alpha_max
        flag = true;
        return;
    end
end

end

%% =========================================================================
% characteristic quadratic interpolation for circular particle
%% =========================================================================

function [ok, val] = characteristic_quadric_value_circle( ...
    phi, P, m, C, r, h, ghostnum, Ny, Nx, start, endx, endy)

ok  = false;
val = NaN;

eps_div = 1e-14;
tol     = 1e-12;

m = m(:) / (norm(m) + eps_div);
P = P(:);
C = C(:);

cand_s    = [];
cand_type = [];
cand_idx  = [];
cand_D    = zeros(2,0);

n_extra_lines = 8;

for sgn = [-1, 1]

    dir = sgn * m;

    [ok_exit, s_exit] = line_exit_circle(P, dir, C, r);

    if ~ok_exit
        continue;
    end

    Xexit = P + s_exit*dir;

    % Vertical scalar grid lines x = x_i.
    if abs(dir(1)) > 1e-14

        x_exit_idx = Xexit(1)/h + ghostnum + 0.5;

        if dir(1) > 0
            i0 = ceil(x_exit_idx - tol);
            step_i = 1;
        else
            i0 = floor(x_exit_idx + tol);
            step_i = -1;
        end

        for q = 0:n_extra_lines

            ii = i0 + q*step_i;

            if ii < 1 || ii > Nx
                continue;
            end

            x_line = (ii - ghostnum - 0.5)*h;
            s = (x_line - P(1)) / dir(1);

            if s < s_exit - 1e-12 || s <= 1e-12
                continue;
            end

            y_hit = P(2) + s*dir(2);
            D = [x_line; y_hit];

            if ~point_inside_active_box(D, h, ghostnum, start, endx, endy)
                continue;
            end

            if norm(D - C) < r - 1e-10
                continue;
            end

            cand_s(end+1)    = s;  %#ok<AGROW>
            cand_type(end+1) = 1;  %#ok<AGROW>
            cand_idx(end+1)  = ii; %#ok<AGROW>
            cand_D(:,end+1)  = D;  %#ok<AGROW>
        end
    end

    % Horizontal scalar grid lines y = y_j.
    if abs(dir(2)) > 1e-14

        y_exit_idx = Xexit(2)/h + ghostnum + 0.5;

        if dir(2) > 0
            j0 = ceil(y_exit_idx - tol);
            step_j = 1;
        else
            j0 = floor(y_exit_idx + tol);
            step_j = -1;
        end

        for q = 0:n_extra_lines

            jj = j0 + q*step_j;

            if jj < 1 || jj > Ny
                continue;
            end

            y_line = (jj - ghostnum - 0.5)*h;
            s = (y_line - P(2)) / dir(2);

            if s < s_exit - 1e-12 || s <= 1e-12
                continue;
            end

            x_hit = P(1) + s*dir(1);
            D = [x_hit; y_line];

            if ~point_inside_active_box(D, h, ghostnum, start, endx, endy)
                continue;
            end

            if norm(D - C) < r - 1e-10
                continue;
            end

            cand_s(end+1)    = s;  %#ok<AGROW>
            cand_type(end+1) = 2;  %#ok<AGROW>
            cand_idx(end+1)  = jj; %#ok<AGROW>
            cand_D(:,end+1)  = D;  %#ok<AGROW>
        end
    end
end

if isempty(cand_s)
    return;
end

[~, order] = sort(cand_s, 'ascend');

for kk = 1:length(order)

    id = order(kk);
    D  = cand_D(:,id);

    if cand_type(id) == 1
        i_line = cand_idx(id);

        [ok_tmp, val_tmp] = quadric_interp_vertical_line_outside_circle( ...
            phi, D, i_line, C, r, h, ghostnum, Ny, Nx);

    else
        j_line = cand_idx(id);

        [ok_tmp, val_tmp] = quadric_interp_horizontal_line_outside_circle( ...
            phi, D, j_line, C, r, h, ghostnum, Ny, Nx);
    end

    if ok_tmp
        ok  = true;
        val = min(max(val_tmp, 0.0), 1.0);
        return;
    end
end

end

function [ok, s_exit] = line_exit_circle(P, dir, C, r)

ok = false;
s_exit = NaN;

Q = P(:) - C(:);
dir = dir(:) / (norm(dir) + 1e-14);

b = dot(Q, dir);
c = dot(Q, Q) - r^2;

disc = b^2 - c;

if disc < 0
    return;
end

root = sqrt(max(disc, 0.0));

s1 = -b - root;
s2 = -b + root;

cands = [s1, s2];
cands = cands(cands > 1e-12);

if isempty(cands)
    return;
end

s_exit = min(cands);
ok = true;

end

function tf = point_inside_active_box(D, h, ghostnum, start, endx, endy)

x = D(1);
y = D(2);

x_min = (start - ghostnum - 0.5)*h;
x_max = (endx  - ghostnum - 0.5)*h;
y_min = (start - ghostnum - 0.5)*h;
y_max = (endy  - ghostnum - 0.5)*h;

tf = (x >= x_min && x <= x_max && y >= y_min && y <= y_max);

end

%% =========================================================================
% quadratic interpolation on vertical scalar grid line
%% =========================================================================

function [ok, val] = quadric_interp_vertical_line_outside_circle( ...
    phi, D, i_line, C, r, h, ghostnum, Ny, Nx)

ok  = false;
val = NaN;

if i_line < 1 || i_line > Nx
    return;
end

x_line = (i_line - ghostnum - 0.5)*h;

eta = D(2)/h + ghostnum + 0.5;

j_near = round(eta);
if j_near >= 1 && j_near <= Ny
    y_near = (j_near - ghostnum - 0.5)*h;

    if abs(D(2) - y_near) < 1e-12
        if is_outside_circle([x_line; y_near], C, r)
            val = phi(j_near, i_line);
            ok  = true;
            return;
        end
    end
end

jA = floor(eta);
jB = jA + 1;

if jA < 1 || jB > Ny
    return;
end

yA = (jA - ghostnum - 0.5)*h;
yB = (jB - ghostnum - 0.5)*h;

A = [x_line; yA];
B = [x_line; yB];

if ~is_outside_circle(A, C, r) || ~is_outside_circle(B, C, r)
    return;
end

cand = [];

jC1 = jA - 1;
if jC1 >= 1
    yC1 = (jC1 - ghostnum - 0.5)*h;
    C1 = [x_line; yC1];

    if is_outside_circle(C1, C, r)
        cand(end+1) = jC1; %#ok<AGROW>
    end
end

jC2 = jB + 1;
if jC2 <= Ny
    yC2 = (jC2 - ghostnum - 0.5)*h;
    C2 = [x_line; yC2];

    if is_outside_circle(C2, C, r)
        cand(end+1) = jC2; %#ok<AGROW>
    end
end

if isempty(cand)
    return;
end

jC = choose_point_nearer_to_circle_boundary_vertical(cand, x_line, C, r, h, ghostnum);
yC = (jC - ghostnum - 0.5)*h;

x_nodes = [yA; yB; yC];
f_nodes = [phi(jA, i_line); phi(jB, i_line); phi(jC, i_line)];

val = lagrange3_eval(x_nodes, f_nodes, D(2));
ok  = true;

end

function jC = choose_point_nearer_to_circle_boundary_vertical(cand, x_line, C, r, h, ghostnum)

jC = cand(1);
best_dist = inf;

for q = 1:length(cand)
    jj = cand(q);
    yq = (jj - ghostnum - 0.5)*h;
    rr = hypot(x_line - C(1), yq - C(2));
    dist_to_boundary = abs(rr - r);

    if dist_to_boundary < best_dist
        best_dist = dist_to_boundary;
        jC = jj;
    end
end

end

%% =========================================================================
% quadratic interpolation on horizontal scalar grid line
%% =========================================================================

function [ok, val] = quadric_interp_horizontal_line_outside_circle( ...
    phi, D, j_line, C, r, h, ghostnum, Ny, Nx)

ok  = false;
val = NaN;

if j_line < 1 || j_line > Ny
    return;
end

y_line = (j_line - ghostnum - 0.5)*h;

xi = D(1)/h + ghostnum + 0.5;

i_near = round(xi);
if i_near >= 1 && i_near <= Nx
    x_near = (i_near - ghostnum - 0.5)*h;

    if abs(D(1) - x_near) < 1e-12
        if is_outside_circle([x_near; y_line], C, r)
            val = phi(j_line, i_near);
            ok  = true;
            return;
        end
    end
end

iA = floor(xi);
iB = iA + 1;

if iA < 1 || iB > Nx
    return;
end

xA = (iA - ghostnum - 0.5)*h;
xB = (iB - ghostnum - 0.5)*h;

A = [xA; y_line];
B = [xB; y_line];

if ~is_outside_circle(A, C, r) || ~is_outside_circle(B, C, r)
    return;
end

cand = [];

iC1 = iA - 1;
if iC1 >= 1
    xC1 = (iC1 - ghostnum - 0.5)*h;
    C1 = [xC1; y_line];

    if is_outside_circle(C1, C, r)
        cand(end+1) = iC1; %#ok<AGROW>
    end
end

iC2 = iB + 1;
if iC2 <= Nx
    xC2 = (iC2 - ghostnum - 0.5)*h;
    C2 = [xC2; y_line];

    if is_outside_circle(C2, C, r)
        cand(end+1) = iC2; %#ok<AGROW>
    end
end

if isempty(cand)
    return;
end

iC = choose_point_nearer_to_circle_boundary_horizontal(cand, y_line, C, r, h, ghostnum);
xC = (iC - ghostnum - 0.5)*h;

x_nodes = [xA; xB; xC];
f_nodes = [phi(j_line, iA); phi(j_line, iB); phi(j_line, iC)];

val = lagrange3_eval(x_nodes, f_nodes, D(1));
ok  = true;

end

function iC = choose_point_nearer_to_circle_boundary_horizontal(cand, y_line, C, r, h, ghostnum)

iC = cand(1);
best_dist = inf;

for q = 1:length(cand)
    ii = cand(q);
    xq = (ii - ghostnum - 0.5)*h;
    rr = hypot(xq - C(1), y_line - C(2));
    dist_to_boundary = abs(rr - r);

    if dist_to_boundary < best_dist
        best_dist = dist_to_boundary;
        iC = ii;
    end
end

end

%% =========================================================================
% small helpers
%% =========================================================================

function tf = is_outside_circle(P, C, r)

tf = (norm(P(:) - C(:)) >= r - 1e-12);

end

function val = lagrange3_eval(x_nodes, f_nodes, xq)

x1 = x_nodes(1);
x2 = x_nodes(2);
x3 = x_nodes(3);

f1 = f_nodes(1);
f2 = f_nodes(2);
f3 = f_nodes(3);

den1 = (x1 - x2) * (x1 - x3);
den2 = (x2 - x1) * (x2 - x3);
den3 = (x3 - x1) * (x3 - x2);

if abs(den1) < 1e-30 || abs(den2) < 1e-30 || abs(den3) < 1e-30
    val = f2;
    return;
end

L1 = ((xq - x2) * (xq - x3)) / den1;
L2 = ((xq - x1) * (xq - x3)) / den2;
L3 = ((xq - x1) * (xq - x2)) / den3;

val = L1*f1 + L2*f2 + L3*f3;
val = min(max(val, 0.0), 1.0);

end

function [ok, D] = first_outside_point_simple(P, m, C, r)

ok = false;
D  = [NaN; NaN];

m = m(:) / (norm(m) + 1e-14);
P = P(:);
C = C(:);

best_s = inf;

for sgn = [-1, 1]

    dir = sgn * m;

    [ok_exit, s_exit] = line_exit_circle(P, dir, C, r);

    if ok_exit && s_exit > 0 && s_exit < best_s
        best_s = s_exit;
        D = P + s_exit * dir;
        ok = true;
    end
end

end