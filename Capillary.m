function [F_cap_x, F_cap_y, T_cap] = Capillary(Fluid, Particle, grid)
%% ========================================================================
% Capillary
%
% Solid-side capillary force / torque calculation.
%
% Modified version for unresolved particle-bubble liquid film:
%
%   1) If grid.enable_alpha_ext_for_capillary == true and Fluid.alpha_ext
%      exists, this function uses Fluid.alpha_ext as the capillary-geometry
%      field.
%
%   2) If Particle.virtual_contact_active(p) == true, alpha=0.5 crossing is
%      searched on the virtual particle boundary:
%
%           r_sample = r_real + grid.virtual_film_thickness
%
%      This allows the particle to feel capillary force after induction /
%      attachment even when the true VOF interface does not touch the real
%      particle boundary because the unresolved film is thinner than grid
%      resolution.
%
%   3) Force direction is computed from the virtual contact-line geometry.
%
%   4) Torque arm is projected back to the real particle surface:
%
%           X_p_force = x_c + r_real*cos(phi_p)
%           Y_p_force = y_c + r_real*sin(phi_p)
%
%      This avoids artificially increasing torque due to the virtual radius.
%
% IMPORTANT:
%   - This function only returns the capillary force/torque on particles.
%   - It does not modify Fluid.alpha or Fluid.alpha_ext.
%   - Fluid density/viscosity should still be based on Fluid.alpha, not
%     Fluid.alpha_ext.
%
% INPUT:
%   Fluid.alpha      : real VOF field
%   Fluid.alpha_ext  : optional capillary-geometry field
%   Particle         : needs x_c, y_c, r, theta/contact_angle/theta_eq
%   grid             : usual grid structure
%
% OUTPUT:
%   F_cap_x, F_cap_y : total capillary force on each particle
%   T_cap            : total capillary torque on each particle
%% ========================================================================

Ny       = grid.Ny;
Nx       = grid.Nx;
h        = grid.h;
ghostnum = grid.ghostnum;

sigma = Fluid.sigma;

Np = length(Particle.r);

F_cap_x = zeros(Np,1);
F_cap_y = zeros(Np,1);
T_cap   = zeros(Np,1);

%% ---------------- global switches ----------------
use_capillary = true;
if isfield(grid, 'enable_capillary')
    use_capillary = logical(grid.enable_capillary);
end

if ~use_capillary || Np == 0
    return;
end

% Whether capillary force should use Fluid.alpha_ext.
use_alpha_ext_for_capillary = false;
if isfield(grid, 'enable_alpha_ext_for_capillary')
    use_alpha_ext_for_capillary = logical(grid.enable_alpha_ext_for_capillary);
end

if use_alpha_ext_for_capillary && isfield(Fluid, 'alpha_ext') && ~isempty(Fluid.alpha_ext)
    alpha = Fluid.alpha_ext;
else
    alpha = Fluid.alpha;
end

alpha = apply_boundary_conditions('vof', grid, alpha);

%% ---------------- parameters ----------------
eps_div = 1e-14;

% Lagrangian resolution along the circular boundary.
N_lag_default = 100;
if isfield(grid, 'capillary_N_lag')
    N_lag_default = grid.capillary_N_lag;
end

% In 2D Cartesian, force is per unit depth => l = 1.
% For axisymmetric cases, use l = 2*pi*r_contact.
is_axisymmetric = false;
if isfield(grid, 'is_axisymmetric')
    is_axisymmetric = logical(grid.is_axisymmetric);
end

% Virtual film thickness.
virtual_film_thickness = 1.5*h;
if isfield(grid, 'virtual_film_thickness')
    virtual_film_thickness = grid.virtual_film_thickness;
end

% If true, only active virtual-contact particles use virtual radius.
% If false, all particles are treated as original real-boundary force.
enable_virtual_capillary = true;
if isfield(grid, 'enable_virtual_capillary')
    enable_virtual_capillary = logical(grid.enable_virtual_capillary);
end

% Interface crossing threshold.
alpha_if = 0.5;
if isfield(grid, 'capillary_alpha_interface')
    alpha_if = grid.capillary_alpha_interface;
end

% Weak filter for numerical tiny crossings.
alpha_min = 1e-3;
alpha_max = 1 - 1e-3;

if isfield(grid, 'capillary_alpha_min')
    alpha_min = grid.capillary_alpha_min;
end
if isfield(grid, 'capillary_alpha_max')
    alpha_max = grid.capillary_alpha_max;
end

%% ========================================================================
% loop over particles
%% ========================================================================
for p = 1:Np

    x_c    = Particle.x_c(p);
    y_c    = Particle.y_c(p);
    r_real = Particle.r(p);

    theta_eq = get_particle_contact_angle_local(Particle, p);   % rad

    % ---------------------------------------------------------------------
    % Decide whether this particle uses real or virtual sampling radius.
    % ---------------------------------------------------------------------
    use_virtual_radius = false;

    if enable_virtual_capillary
        if isfield(Particle, 'virtual_contact_active') && ...
                numel(Particle.virtual_contact_active) >= p && ...
                Particle.virtual_contact_active(p)
            use_virtual_radius = true;
        end
    end

    if use_virtual_radius
        r_sample = r_real + virtual_film_thickness;
    else
        r_sample = r_real;
    end

    % ---------------------------------------------------------------------
    % Discretize the sampling circular boundary.
    % For active virtual contact, this is the virtual particle boundary.
    % Otherwise, this is the real particle boundary.
    % ---------------------------------------------------------------------
    N_lag = max(N_lag_default, ceil(2*pi*r_sample/(0.5*h)));

    phi = linspace(0, 2*pi, N_lag+1);
    phi = phi(1:end-1);

    alpha_bnd = zeros(1, N_lag);

    % ---------------------------------------------------------------------
    % Step 1: sample alpha on selected particle boundary
    % ---------------------------------------------------------------------
    for k = 1:N_lag

        xb = x_c + r_sample * cos(phi(k));
        yb = y_c + r_sample * sin(phi(k));

        alpha_bnd(k) = interp_cell_scalar_bilinear_local( ...
            alpha, xb, yb, h, ghostnum, Ny, Nx);
    end

    % ---------------------------------------------------------------------
    % Step 2: find alpha = alpha_if crossings on selected boundary
    % ---------------------------------------------------------------------
    for k = 1:N_lag

        kp = k + 1;
        if kp > N_lag
            kp = 1;
        end

        gamma_A = alpha_bnd(k);
        gamma_B = alpha_bnd(kp);

        near_contact = ...
            ((gamma_A > alpha_min && gamma_A < alpha_max) || ...
             (gamma_B > alpha_min && gamma_B < alpha_max));

        if ~near_contact
            continue;
        end

        % Need a genuine crossing.
        if (gamma_A - alpha_if) * (gamma_B - alpha_if) > 0
            continue;
        end

        if abs(gamma_A - gamma_B) < 1e-10
            continue;
        end

        phi_A = phi(k);
        phi_B = phi(kp);

        % Handle 2pi -> 0 branch cut.
        if kp == 1 && phi_B < phi_A
            phi_B = phi_B + 2*pi;
        end

        % -----------------------------------------------------------------
        % Step 3: angular interpolation to find contact location.
        % -----------------------------------------------------------------
        lambda = (alpha_if - gamma_A) / (gamma_B - gamma_A + eps_div);
        lambda = min(max(lambda, 0.0), 1.0);

        dphi  = lambda * (phi_B - phi_A);
        phi_p = phi_A + dphi;

        % Normalize for trig use.
        phi_p_mod = mod(phi_p, 2*pi);

        % -----------------------------------------------------------------
        % Step 4: contact point on the sampling boundary.
        %
        % X_p_sample:
        %   position used for geometric direction construction.
        %
        % X_p_force:
        %   real-solid projected position used for torque arm.
        % -----------------------------------------------------------------
        X_p_sample = x_c + r_sample * cos(phi_p_mod);
        Y_p_sample = y_c + r_sample * sin(phi_p_mod);

        X_p_force = x_c + r_real * cos(phi_p_mod);
        Y_p_force = y_c + r_real * sin(phi_p_mod);

        % Point A on sampling boundary.
        X_A = x_c + r_sample * cos(phi_A);
        Y_A = y_c + r_sample * sin(phi_A);

        %#ok<NASGU>
        dummy_XA = X_A;
        dummy_YA = Y_A;

        % -----------------------------------------------------------------
        % Step 5: construct t_A and rotate to get t_cl.
        %
        % For circular boundary:
        %   outward normal at A:
        %       n_A = [cos(phi_A); sin(phi_A)]
        %
        %   CCW tangent at A:
        %       t_wA = [-sin(phi_A); cos(phi_A)]
        %
        % Boundary alpha variation determines the liquid-side tangential
        % direction.
        % -----------------------------------------------------------------
        n_A  = [cos(phi_A); sin(phi_A)];
        t_wA = [-sin(phi_A); cos(phi_A)];

        % If gamma increases from A -> B, take +t_wA as liquid-side direction.
        % Otherwise take -t_wA.
        if gamma_B > gamma_A
            t_liq_A = t_wA;
        else
            t_liq_A = -t_wA;
        end

        % Two interface-tangent candidates obtained from contact angle.
        % Keep the one pointing away from the solid.
        tA_cand1 = rotate_vec_local(t_liq_A, +theta_eq);
        tA_cand2 = rotate_vec_local(t_liq_A, -theta_eq);

        s1 = dot(tA_cand1, n_A);
        s2 = dot(tA_cand2, n_A);

        if s1 >= s2
            t_A = tA_cand1;
        else
            t_A = tA_cand2;
        end

        t_A = t_A / (norm(t_A) + eps_div);

        % -----------------------------------------------------------------
        % Step 6: rotate t_A by dphi to obtain t_cl.
        % -----------------------------------------------------------------
        t_cl = rotate_vec_local(t_A, dphi);
        t_cl = t_cl / (norm(t_cl) + eps_div);

        % -----------------------------------------------------------------
        % Step 7: contact-line length and Young force.
        %
        % For axisymmetric cases, use the real projected contact position
        % for the ring length to avoid artificial virtual-radius torque/force
        % amplification.
        % -----------------------------------------------------------------
        if is_axisymmetric
            % Assuming x is radial coordinate.
            l_cl = 2*pi*abs(X_p_force);
        else
            % 2D Cartesian per unit depth.
            l_cl = 1.0;
        end

        F_sigma = sigma * l_cl * t_cl;

        F_cap_x(p) = F_cap_x(p) + F_sigma(1);
        F_cap_y(p) = F_cap_y(p) + F_sigma(2);

        % Torque about real particle center using real projected contact arm.
        rx = X_p_force - x_c;
        ry = Y_p_force - y_c;

        T_cap(p) = T_cap(p) + (rx * F_sigma(2) - ry * F_sigma(1));

        % Keep sample point variables used, avoiding accidental removal in
        % future debugging.
        
        dummy_sample_x = X_p_sample;
        dummy_sample_y = Y_p_sample;
    end
end

end

%% ========================================================================
% helper: rotate 2D vector
%% ========================================================================
function vout = rotate_vec_local(vin, ang)

c = cos(ang);
s = sin(ang);

R = [c, -s; ...
     s,  c];

vout = R * vin(:);

end

%% ========================================================================
% helper: bilinear interpolation of cell-centered scalar
%% ========================================================================
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

%% ========================================================================
% helper: particle contact angle
%% ========================================================================
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