function [u_pred, v_pred] = prediction(Fluid, grid)
%% ========================================================================
% prediction: height-function CSF + unified face-centered 1/rho
%
% Main features:
%   1) Uses compute_face_properties() for invrho_x/invrho_y.
%   2) Pressure gradient uses invrho_x/invrho_y.
%   3) Surface tension uses the same invrho_x/invrho_y.
%   4) Curvature is computed by a height-function method where available.
%   5) If grid.fix_bubble == true, fluid CSF is disabled.
%
% Projection convention:
%   This prediction step uses old pressure p^n.
%   Poisson solves pressure correction phi.
%   p^{n+1} = p^n + phi.
%
% Main modification:
%   In CSF, the default alpha gradient uses smoothed alpha_cv instead of
%   raw alpha. This improves consistency between curvature and interface
%   force localization and usually reduces parasitic currents.
%% ========================================================================

start = grid.start;
endy  = grid.endy;
endx  = grid.endx;
Ny    = grid.Ny;
Nx    = grid.Nx;
h     = grid.h;
dt    = grid.dt;

u_in = Fluid.u;
v_in = Fluid.v;
p_in = Fluid.p;

% -------------------------------------------------------------------------
% Real alpha and CSF alpha
%
% alpha_real:
%   Real VOF field. Used for density, viscosity, mass conservation, etc.
%
% alpha_csf:
%   Auxiliary capillary-geometry field. If Fluid.alpha_ext exists and
%   grid.enable_alpha_ext_for_csf == true, CSF uses alpha_ext.
%   Otherwise CSF falls back to real alpha.
% -------------------------------------------------------------------------
alpha_real = Fluid.alpha;
alpha_real = apply_boundary_conditions('vof', grid, alpha_real);

use_alpha_ext_for_csf = false;
if isfield(grid, 'enable_alpha_ext_for_csf')
    use_alpha_ext_for_csf = logical(grid.enable_alpha_ext_for_csf);
end

if use_alpha_ext_for_csf && isfield(Fluid, 'alpha_ext') && ~isempty(Fluid.alpha_ext)
    alpha_csf = Fluid.alpha_ext;
else
    alpha_csf = alpha_real;
end

alpha_csf = apply_boundary_conditions('vof', grid, alpha_csf);

sigma = Fluid.sigma;

% Fluid gravity removed in your current framework.
if isfield(grid, 'g')
    g = grid.g;
else
    g = 9.81;
end

%% ---------------- unified face properties ----------------
prop = compute_face_properties(Fluid, grid);

invrho_x = prop.invrho_x;
invrho_y = prop.invrho_y;
mu_x     = prop.mu_x;
mu_y     = prop.mu_y;

%% ---------------- advection ----------------
adux = zeros(Ny+1, Nx+1);
aduy = zeros(Ny+1, Nx+1);
advx = zeros(Ny+1, Nx+1);
advy = zeros(Ny+1, Nx+1);

for j = start:endy+1
    for i = start:endx+1

        %% du/dx
        if u_in(j,i) > 0
            adux(j,i) = u_in(j,i) * ...
                (3*u_in(j,i) - 4*u_in(j,i-1) + u_in(j,i-2)) / (2*h);
        else
            adux(j,i) = u_in(j,i) * ...
                (-3*u_in(j,i) + 4*u_in(j,i+1) - u_in(j,i+2)) / (2*h);
        end

        %% du/dy
        v_at_u = 0.25 * ...
            (v_in(j,i-1) + v_in(j,i) + v_in(j+1,i-1) + v_in(j+1,i));

        if v_at_u > 0
            aduy(j,i) = v_at_u * ...
                (3*u_in(j,i) - 4*u_in(j-1,i) + u_in(j-2,i)) / (2*h);
        else
            aduy(j,i) = v_at_u * ...
                (-3*u_in(j,i) + 4*u_in(j+1,i) - u_in(j+2,i)) / (2*h);
        end
    end
end

adu = adux + aduy;

for j = start:endy+1
    for i = start:endx+1

        %% dv/dx
        u_at_v = 0.25 * ...
            (u_in(j-1,i) + u_in(j-1,i+1) + u_in(j,i) + u_in(j,i+1));

        if u_at_v > 0
            advx(j,i) = u_at_v * ...
                (3*v_in(j,i) - 4*v_in(j,i-1) + v_in(j,i-2)) / (2*h);
        else
            advx(j,i) = u_at_v * ...
                (-3*v_in(j,i) + 4*v_in(j,i+1) - v_in(j,i+2)) / (2*h);
        end

        %% dv/dy
        if v_in(j,i) > 0
            advy(j,i) = v_in(j,i) * ...
                (3*v_in(j,i) - 4*v_in(j-1,i) + v_in(j-2,i)) / (2*h);
        else
            advy(j,i) = v_in(j,i) * ...
                (-3*v_in(j,i) + 4*v_in(j+1,i) - v_in(j+2,i)) / (2*h);
        end
    end
end

adv = advx + advy;

%% ---------------- diffusion ----------------
diffu = zeros(Ny+1, Nx+1);
diffv = zeros(Ny+1, Nx+1);

for j = start:endy+1
    for i = start:endx+1

        lap_u = ...
            (u_in(j+1,i) - 2*u_in(j,i) + u_in(j-1,i)) / h^2 + ...
            (u_in(j,i+1) - 2*u_in(j,i) + u_in(j,i-1)) / h^2;

        diffu(j,i) = mu_x(j,i) * lap_u * invrho_x(j,i);

        lap_v = ...
            (v_in(j+1,i) - 2*v_in(j,i) + v_in(j-1,i)) / h^2 + ...
            (v_in(j,i+1) - 2*v_in(j,i) + v_in(j,i-1)) / h^2;

        diffv(j,i) = mu_y(j,i) * lap_v * invrho_y(j,i);
    end
end

%% ---------------- height-function CSF ----------------
[surfu, surfv] = surface_tension_height_function( ...
    alpha_csf, sigma, h, Ny, Nx, start, endy, endx, ...
    invrho_x, invrho_y, grid);



% In fixed-bubble tests, do not let fluid CSF generate parasitic currents.
% Particle capillary force from Capillary() is not affected.
if isfield(grid, 'fix_bubble') && grid.fix_bubble
    surfu(:) = 0;
    surfv(:) = 0;
end

%% ---------------- pressure gradient using same invrho faces ----------------
Pu = zeros(Ny+1, Nx+1);
Pv = zeros(Ny+1, Nx+1);

for j = start:endy+1
    for i = start:endx+1
        Pu(j,i) = invrho_x(j,i) * (p_in(j,i) - p_in(j,i-1)) / h;
        Pv(j,i) = invrho_y(j,i) * (p_in(j,i) - p_in(j-1,i)) / h;
    end
end

%% ---------------- update ----------------
u_pred = u_in + (-adu + diffu + surfu - Pu) * dt;
v_pred = v_in + (-adv + diffv + surfv - Pv - g) * dt;

[u_pred, v_pred] = apply_boundary_conditions('velocity', grid, u_pred, v_pred);

end


function [surfu, surfv, kappa, hf_ok, kappa_dn] = surface_tension_height_function( ...
    alpha, sigma, h, Ny, Nx, start, endy, endx, ...
    invrho_x, invrho_y, grid)

%% ========================================================================
% Height-function based CSF without 45-degree special repair
% and without global curvature sign forcing.
%
% Main changes relative to the previous version:
%   1) Removed all 45-degree patch repair logic.
%   2) Removed global curvature sign forcing.
%   3) Height-function curvature keeps local sign from kappa_dn = -div(n).
%      Therefore dimple regions are allowed to have opposite signed curvature.
%   4) Retains:
%        - smoothed alpha for orientation/interface detection;
%        - raw alpha for height-function volume sums;
%        - neighbor fill for failed HF points;
%        - weak interface-only smoothing;
%        - well-balanced capillary-pressure surface-tension force:
%              q = sigma * kappa * alpha
%              F = grad(q)
%
% Outputs:
%   surfu     surface-tension acceleration at u-faces
%   surfv     surface-tension acceleration at v-faces
%   kappa     final cell-centered signed curvature
%   hf_ok     height-function valid mask
%   kappa_dn  fallback curvature from -div(n)
%% ========================================================================

eps_if = 1e-3;

% Curvature clipping. Keep this configurable.
max_k = 1.0 / h;
if isfield(grid, 'csf_kappa_clip_factor')
    max_k = grid.csf_kappa_clip_factor / h;
end

%% ========================================================================
% Controls
%% ========================================================================

% Fill failed HF cells from neighboring successful HF cells.
enable_hf_neighbor_fill = true;
if isfield(grid, 'csf_hf_neighbor_fill')
    enable_hf_neighbor_fill = grid.csf_hf_neighbor_fill;
end

% Weak curvature smoothing after local repairs.
enable_kappa_smoothing = false;
if isfield(grid, 'csf_kappa_smoothing')
    enable_kappa_smoothing = grid.csf_kappa_smoothing;
end

num_smooth_iter = 1;
if isfield(grid, 'csf_kappa_smooth_iter')
    num_smooth_iter = grid.csf_kappa_smooth_iter;
end

%% ========================================================================
% 1. Smooth alpha for orientation and interface detection only
%% ========================================================================
alpha_cv = alpha;

for j = start:endy
    for i = start:endx
        alpha_cv(j,i) = ( ...
            1*alpha(j-1,i-1) + 2*alpha(j-1,i) + 1*alpha(j-1,i+1) + ...
            2*alpha(j  ,i-1) + 4*alpha(j  ,i) + 2*alpha(j  ,i+1) + ...
            1*alpha(j+1,i-1) + 2*alpha(j+1,i) + 1*alpha(j+1,i+1) ) / 16;
    end
end

alpha_cv = apply_boundary_conditions('vof', grid, alpha_cv);

%% Interface mask
interface_mask = false(Ny, Nx);

for j = start:endy
    for i = start:endx
        if alpha_cv(j,i) > eps_if && alpha_cv(j,i) < 1 - eps_if
            interface_mask(j,i) = true;
        end
    end
end

%% ========================================================================
% 2. Fallback curvature from -div(n)
%
% Here n = grad(alpha_cv) / |grad(alpha_cv)|.
% kappa_dn is signed and is used as the local sign reference for HF curvature.
%% ========================================================================
nx = zeros(Ny, Nx);
ny = zeros(Ny, Nx);

for j = start:endy
    for i = start:endx
        ax = (alpha_cv(j,i+1) - alpha_cv(j,i-1)) / (2*h);
        ay = (alpha_cv(j+1,i) - alpha_cv(j-1,i)) / (2*h);

        gmag = sqrt(ax^2 + ay^2) + 1e-20;

        nx(j,i) = ax / gmag;
        ny(j,i) = ay / gmag;
    end
end

kappa_dn = zeros(Ny, Nx);

for j = start:endy
    for i = start:endx

        kval = - ( ...
            (nx(j,i+1) - nx(j,i-1)) / (2*h) + ...
            (ny(j+1,i) - ny(j-1,i)) / (2*h) );

        if abs(kval) > max_k
            kval = sign(kval) * max_k;
        end

        kappa_dn(j,i) = kval;
    end
end

kappa_dn = apply_boundary_conditions('pressure', grid, kappa_dn);

%% ========================================================================
% 3. Height-function curvature
%
% No global sign forcing.
% The HF formula is used for the magnitude, and its sign is taken from the
% local signed fallback curvature kappa_dn. If local kappa_dn is too small,
% use nearby signed fallback values.
%% ========================================================================
kappa = kappa_dn;
kappa_hf = zeros(Ny, Nx);
hf_ok = false(Ny, Nx);

r = 3;

for j = start+r:endy-r
    for i = start+r:endx-r

        if ~interface_mask(j,i)
            continue;
        end

        if abs(ny(j,i)) >= abs(nx(j,i))
            [ok, kval_abs] = hf_curvature_vertical(alpha, i, j, h, r);
        else
            [ok, kval_abs] = hf_curvature_horizontal(alpha, i, j, h, r);
        end

        if ok
            local_sign = get_local_kappa_sign(kappa_dn, interface_mask, i, j, start, endx, endy);

            kval = local_sign * abs(kval_abs);

            if abs(kval) > max_k
                kval = sign(kval) * max_k;
            end

            kappa_hf(j,i) = kval;
            kappa(j,i) = kval;
            hf_ok(j,i) = true;
        end
    end
end

%% ========================================================================
% 4. Replace isolated HF-failed interface points using neighboring HF points
%% ========================================================================
if enable_hf_neighbor_fill

    kappa_filled = kappa;
    hf_filled = hf_ok;

    % First pass: fill failed interface cells with neighboring HF curvature.
    for j = start+1:endy-1
        for i = start+1:endx-1

            if interface_mask(j,i) && ~hf_ok(j,i)

                vals = [];
                weights = [];

                for jj = j-1:j+1
                    for ii = i-1:i+1
                        if hf_ok(jj,ii)

                            dist2 = (jj-j)^2 + (ii-i)^2;
                            w_dist = 1 / (dist2 + 1);

                            w_if = alpha_cv(jj,ii) * (1 - alpha_cv(jj,ii));
                            w = w_dist * (w_if + 1e-12);

                            vals(end+1) = kappa_hf(jj,ii); %#ok<AGROW>
                            weights(end+1) = w; %#ok<AGROW>
                        end
                    end
                end

                if ~isempty(vals)
                    kappa_filled(j,i) = sum(weights .* vals) / (sum(weights) + 1e-30);
                    hf_filled(j,i) = true;
                end
            end
        end
    end

    % Second pass: if still failed, use a wider 5x5 neighborhood.
    for j = start+2:endy-2
        for i = start+2:endx-2

            if interface_mask(j,i) && ~hf_filled(j,i)

                vals = [];
                weights = [];

                for jj = j-2:j+2
                    for ii = i-2:i+2
                        if hf_ok(jj,ii)

                            dist2 = (jj-j)^2 + (ii-i)^2;
                            w_dist = 1 / (dist2 + 1);

                            w_if = alpha_cv(jj,ii) * (1 - alpha_cv(jj,ii));
                            w = w_dist * (w_if + 1e-12);

                            vals(end+1) = kappa_hf(jj,ii); %#ok<AGROW>
                            weights(end+1) = w; %#ok<AGROW>
                        end
                    end
                end

                if ~isempty(vals)
                    kappa_filled(j,i) = sum(weights .* vals) / (sum(weights) + 1e-30);
                    hf_filled(j,i) = true;
                end
            end
        end
    end

    kappa = kappa_filled;
    hf_ok = hf_filled;
end

%% ========================================================================
% 5. For remaining failed interface cells, use smoothed fallback curvature
%
% No global sign forcing here. The fallback curvature is already signed.
%% ========================================================================
kappa_tmp = kappa;

for j = start+1:endy-1
    for i = start+1:endx-1

        if interface_mask(j,i) && ~hf_ok(j,i)

            vals = [];
            weights = [];

            for jj = j-1:j+1
                for ii = i-1:i+1
                    if interface_mask(jj,ii)

                        w_if = alpha_cv(jj,ii) * (1 - alpha_cv(jj,ii));
                        vals(end+1) = kappa_dn(jj,ii); %#ok<AGROW>
                        weights(end+1) = w_if + 1e-12; %#ok<AGROW>
                    end
                end
            end

            if ~isempty(vals)
                kval = sum(weights .* vals) / (sum(weights) + 1e-30);
            else
                kval = kappa_dn(j,i);
            end

            if abs(kval) > max_k
                kval = sign(kval) * max_k;
            end

            kappa_tmp(j,i) = kval;
        end
    end
end

kappa = kappa_tmp;
kappa = apply_boundary_conditions('pressure', grid, kappa);

%% ========================================================================
% 6. Interface-only weak curvature smoothing
%
% No global sign forcing. This allows real signed-curvature reversal in
% dimple regions.
%% ========================================================================
if enable_kappa_smoothing

    for iter_s = 1:num_smooth_iter

        kappa_s = kappa;

        for j = start+1:endy-1
            for i = start+1:endx-1

                if interface_mask(j,i)

                    vals = [];
                    weights = [];

                    for jj = j-1:j+1
                        for ii = i-1:i+1

                            if interface_mask(jj,ii)

                                % Center cell has larger weight.
                                if jj == j && ii == i
                                    w_center = 4;
                                elseif jj == j || ii == i
                                    w_center = 2;
                                else
                                    w_center = 1;
                                end

                                w_if = alpha_cv(jj,ii) * (1 - alpha_cv(jj,ii));
                                w = w_center * (w_if + 1e-12);

                                vals(end+1) = kappa(jj,ii); %#ok<AGROW>
                                weights(end+1) = w; %#ok<AGROW>
                            end
                        end
                    end

                    if ~isempty(vals)
                        kval = sum(weights .* vals) / (sum(weights) + 1e-30);

                        % Do not over-smooth: blend with original value.
                        kappa_s(j,i) = 0.5 * kappa(j,i) + 0.5 * kval;

                        if abs(kappa_s(j,i)) > max_k
                            kappa_s(j,i) = sign(kappa_s(j,i)) * max_k;
                        end
                    end
                end
            end
        end

        kappa = kappa_s;
        kappa = apply_boundary_conditions('pressure', grid, kappa);
    end
end

%% ========================================================================
% 7. Well-balanced capillary-pressure surface tension
%
% Instead of:
%     F = sigma * kappa * grad(alpha)
%
% use:
%     q = sigma * kappa * alpha
%     F = grad(q)
%
% The face gradient is discretized using the same stencil as the pressure
% gradient in prediction:
%
%     Pu = invrho_x * (p(i)-p(i-1))/h
%     Pv = invrho_y * (p(j)-p(j-1))/h
%% ========================================================================

surfu = zeros(Ny+1, Nx+1);
surfv = zeros(Ny+1, Nx+1);

% Use smoothed alpha for the capillary pressure potential.
alpha_force = alpha_cv;

% Only keep q inside the interface band and immediate neighboring cells.
q = zeros(Ny, Nx);

for j = start:endy
    for i = start:endx

        local_interface = false;

        for jj = max(start,j-1):min(endy,j+1)
            for ii = max(start,i-1):min(endx,i+1)
                if alpha_cv(jj,ii) > eps_if && alpha_cv(jj,ii) < 1 - eps_if
                    local_interface = true;
                    break;
                end
            end

            if local_interface
                break;
            end
        end

        if local_interface
            q(j,i) = sigma * kappa(j,i) * alpha_force(j,i);
        else
            q(j,i) = 0;
        end
    end
end

q = apply_boundary_conditions('pressure', grid, q);

% u-face capillary acceleration
for j = start:endy+1
    for i = start:endx+1

        qL = q(j, i-1);
        qR = q(j, i);

        Fux = (qR - qL) / h;

        surfu(j,i) = invrho_x(j,i) * Fux;
    end
end

% v-face capillary acceleration
for j = start:endy+1
    for i = start:endx+1

        qB = q(j-1, i);
        qT = q(j,   i);

        Fvy = (qT - qB) / h;

        surfv(j,i) = invrho_y(j,i) * Fvy;
    end
end

end

%% ========================================================================
% Helper: local curvature sign
%% ========================================================================
function sgn = get_local_kappa_sign(kappa_dn, interface_mask, i, j, start, endx, endy)

eps_k = 1e-12;

if abs(kappa_dn(j,i)) > eps_k && isfinite(kappa_dn(j,i))
    sgn = sign(kappa_dn(j,i));
    return;
end

vals = [];

j1 = max(start, j-1);
j2 = min(endy,  j+1);
i1 = max(start, i-1);
i2 = min(endx,  i+1);

for jj = j1:j2
    for ii = i1:i2
        if interface_mask(jj,ii) && isfinite(kappa_dn(jj,ii)) && abs(kappa_dn(jj,ii)) > eps_k
            vals(end+1) = kappa_dn(jj,ii); %#ok<AGROW>
        end
    end
end

if isempty(vals)
    sgn = 1;
else
    m = mean(vals);
    if abs(m) > eps_k
        sgn = sign(m);
    else
        s = sum(sign(vals));
        if s == 0
            sgn = sign(vals(1));
        else
            sgn = sign(s);
        end
    end
end

if sgn == 0
    sgn = 1;
end

end

%% ========================================================================
% Helper: height-function curvature, horizontal orientation
%% ========================================================================
function [ok, kappa_abs] = hf_curvature_horizontal(alpha, i, j, h, r)
% Interface locally represented as x = H(y).
%
% This function returns curvature magnitude.
% The sign is assigned outside using local kappa_dn.

ok = false;
kappa_abs = 0;

rows = [j-1, j, j+1];
H = zeros(1,3);

for c = 1:3
    jj = rows(c);
    row = alpha(jj, i-r:i+r);

    if ~(min(row) < 0.5 && max(row) > 0.5)
        return;
    end

    H(c) = sum(row) * h;
end

Hy  = (H(3) - H(1)) / (2*h);
Hyy = (H(3) - 2*H(2) + H(1)) / h^2;

kappa_abs = abs(Hyy / ((1 + Hy^2)^(3/2) + 1e-20));

if isfinite(kappa_abs)
    ok = true;
end

end

%% ========================================================================
% Helper: height-function curvature, vertical orientation
%% ========================================================================
function [ok, kappa_abs] = hf_curvature_vertical(alpha, i, j, h, r)
% Interface locally represented as y = H(x).
%
% This function returns curvature magnitude.
% The sign is assigned outside using local kappa_dn.

ok = false;
kappa_abs = 0;

cols = [i-1, i, i+1];
H = zeros(1,3);

for c = 1:3
    ii = cols(c);
    col = alpha(j-r:j+r, ii);

    if ~(min(col) < 0.5 && max(col) > 0.5)
        return;
    end

    H(c) = sum(col) * h;
end

Hx  = (H(3) - H(1)) / (2*h);
Hxx = (H(3) - 2*H(2) + H(1)) / h^2;

kappa_abs = abs(Hxx / ((1 + Hx^2)^(3/2) + 1e-20));

if isfinite(kappa_abs)
    ok = true;
end

end