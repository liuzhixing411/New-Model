clc;
clear;
close all;

%% ========================================================================
%  Single Bubble + Single Particle Virtual Contact Diagnostic Test
%
%  Purpose:
%    1) Fast small-scale test for one bubble + one particle.
%    2) Test particle-bubble state machine:
%          0 free, 1 collision, 2 induction, 3 attachment
%    3) Test separated alpha fields:
%          Fluid.alpha     : real conservative VOF field
%          Fluid.alpha_ext : capillary-geometry extension field
%    4) Test update_extension only after attachment.
%    5) Visualize speed magnitude with contourf and sparse streamlines.
%
%  Required functions:
%    update_particle_bubble_state
%    update_advection
%    update_extension
%    prediction
%    Interaction
%    Capillary
%    update_particle
%    update_fluid
%    apply_boundary_conditions
%    Poisson
%    compute_face_properties
%
%  alpha convention:
%    alpha = 1 : liquid
%    alpha = 0 : gas bubble
%% ========================================================================

%% ---------------- User controls ----------------
t_end = 0.030;
plot_every_nsteps = 4;
diag_every_nsteps = 8;

speed_caxis_max = 0.035;

% For plotting only: do not over-interpolate because this is a fast test.
viz_up = 2;

%% ========================================================================
% 1. Domain and grid
%% ========================================================================
W = 1.50e-3;
L = 1.50e-3;

h = 2.5e-5;
ghostnum = 3;

N = round(W / h);
M = round(L / h);

start = ghostnum + 1;
endy  = M + ghostnum;
endx  = N + ghostnum;

Ny = 2*ghostnum + M;
Nx = 2*ghostnum + N;

Wall.up    = L;
Wall.down  = 0.0;
Wall.left  = 0.0;
Wall.right = W;

% Active cell-center coordinates.
x_vec = ((start:endx) - ghostnum - 0.5) * h;
y_vec = ((start:endy) - ghostnum - 0.5) * h;
[XX, YY] = meshgrid(x_vec, y_vec);

grid.M = M;
grid.N = N;
grid.ghostnum = ghostnum;
grid.start = start;
grid.endy = endy;
grid.endx = endx;
grid.Ny = Ny;
grid.Nx = Nx;
grid.h = h;

% Initial dt. The loop will update it.
grid.dt = 1.0e-5;

% Boundary conditions.
grid.bc.left  = 'free-slip';
grid.bc.right = 'free-slip';
grid.bc.down  = 'free-slip';
grid.bc.up    = 'free-slip';

% prediction.m uses grid.g for fluid body force.
% update_particle.m still has its own particle gravity internally.
grid.g = 0.0;

%% ========================================================================
% 2. Poisson controls
%% ========================================================================
grid.poisson_tol = 1e-4;
grid.poisson_max_iter = 25000;
grid.poisson_omega = 1.35;
grid.poisson_verbose = false;

%% ========================================================================
% 3. Virtual contact / induction controls
%% ========================================================================
grid.enable_particle_bubble_state = true;

% Keep collision state available, but you can close it.
grid.enable_collision_state = true;

% If false, h_p_model <= 0 directly enters attachment.
grid.enable_induction_state = true;

% update_extension is active only when Particle.virtual_contact_active is true.
grid.enable_virtual_extension = true;

% CSF uses Fluid.alpha_ext only if this is true.
% For safer first test, set false; for full coupling, set true.
grid.enable_alpha_ext_for_csf = true;

% Particle capillary force uses Fluid.alpha_ext and virtual radius if true.
grid.enable_alpha_ext_for_capillary = true;
grid.enable_virtual_capillary = true;

grid.enable_capillary = true;

% Virtual unresolved liquid-film thickness.
grid.virtual_film_thickness = 1.5 * h;

% Literature-inspired collision threshold h_c = 0.2 r_p.
grid.collision_hc_factor = 0.2;

% Hysteresis release distance.
grid.virtual_release_thickness = 2.5 * h;

% State detector sampling.
grid.pb_N_lag = 180;
grid.pb_search_ds_factor = 0.35;
grid.pb_search_max_factor_h = 7.0;

% Induction time model:
%   t_a = K * db^2.38 * dp^1.59 / (db+dp)^1.59
% K is unit-dependent; here t_a is clipped for a quick numerical test.
grid.induction_K = 1.0;
grid.induction_t_min = 2 * grid.dt;
grid.induction_t_max = 18 * grid.dt;

% Local bubble radius estimate.
grid.local_Rb_min = 2.0 * h;
grid.local_Rb_max_factor = 1.0;
grid.local_Rb_window = 4.0 * h;

%% ========================================================================
% 4. Capillary / extension / CSF controls
%% ========================================================================
grid.alpha_represents_liquid = true;

grid.capillary_alpha_min = 1.0e-3;
grid.capillary_alpha_max = 1.0 - 1.0e-3;
grid.capillary_probe_dist = 0.5 * h;
grid.capillary_N_lag = 180;

% Real solid inner extension depth.
grid.extension_distance = 2.0 * h;

% Tangential extension width near contact line.
grid.extension_tangent_halfwidth = 3.0 * h;

% For virtual contact, this is usually too strict.
grid.extension_require_interfacial_neighbor = false;

% To make this small test robust near boundaries / coarse grid.
grid.extension_allow_bilinear_fallback = true;

grid.clean_deep_solid_alpha = false;

% CSF controls.
grid.csf_kappa_clip_factor = 1.0;
grid.csf_hf_neighbor_fill = true;
grid.csf_kappa_smoothing = false;
grid.csf_kappa_smooth_iter = 1;

%% ========================================================================
% 5. Fluid initialization
%% ========================================================================
Fluid.rhol  = 1000.0;
Fluid.rhog  = 1.3;
Fluid.sigma = 0.030;
Fluid.mul   = 1.0e-3;
Fluid.mug   = 1.8e-5;

% Staggered velocities: Ny+1 by Nx+1.
Fluid.u = zeros(Ny+1, Nx+1);
Fluid.v = zeros(Ny+1, Nx+1);

% Pressure and alpha are cell-centered: Ny by Nx.
Fluid.p = zeros(Ny, Nx);
Fluid.alpha = ones(Ny, Nx);

%% ---------------- Single bubble ----------------
bub.x = 0.50 * W;
bub.y = 0.72e-3;
bub.r = 2.60e-4;

bubble_area0 = pi * bub.r^2;

% Subcell-sampled initial bubble.
Fluid.alpha = embed_bubbles_union(Fluid.alpha, bub, grid, 24);
Fluid.alpha = fill_vof_ghost_neumann(Fluid.alpha, grid);
Fluid.alpha = apply_boundary_conditions('vof', grid, Fluid.alpha);

% Initial pressure with a simple Laplace jump inside gas bubble.
Fluid.p(:,:) = 0.0;
for j = start:endy
    for i = start:endx
        x = (i - ghostnum - 0.5) * h;
        y = (j - ghostnum - 0.5) * h;

        if hypot(x - bub.x, y - bub.y) <= bub.r
            Fluid.p(j,i) = Fluid.sigma / bub.r;
        end
    end
end
Fluid.p = apply_boundary_conditions('pressure', grid, Fluid.p);

% Auxiliary alpha field.
Fluid.alpha_ext = Fluid.alpha;

[Fluid.u, Fluid.v] = apply_boundary_conditions('velocity', grid, Fluid.u, Fluid.v);

%% ========================================================================
% 6. Particle initialization
%% ========================================================================
Np = 1;

Rp = 5.0e-5;

Particle.x_c = zeros(1,Np);
Particle.y_c = zeros(1,Np);
Particle.r   = zeros(1,Np);
Particle.rho = zeros(1,Np);

% Place particle above bubble, slightly off-center.
initial_gap = 2 * h;

Particle.x_c(1) = bub.x + 0.12 * bub.r;
Particle.y_c(1) = bub.y + bub.r + Rp + initial_gap;
Particle.r(1)   = Rp;

% Light-ish particle for quick diagnostic motion.
Particle.rho(1) = 1500.0;

Particle.u = zeros(1,Np);
Particle.v = -1.5e-2 * ones(1,Np);
Particle.omega = zeros(1,Np);

Particle.V = pi .* Particle.r.^2;
Particle.m = Particle.rho .* Particle.V;
Particle.I = 0.5 .* Particle.r.^2 .* Particle.m;

% Contact angle.
Particle.theta = 110 * ones(1,Np);

% State-machine fields.
Particle.pb_state = zeros(Np,1);
Particle.induction_timer = zeros(Np,1);
Particle.ta_induction = zeros(Np,1);
Particle.Rb_induction = zeros(Np,1);
Particle.virtual_contact_active = false(Np,1);

Particle.pb_contact_phi = nan(Np,1);
Particle.pb_contact_x = nan(Np,1);
Particle.pb_contact_y = nan(Np,1);
Particle.pb_hp_geo = inf(Np,1);
Particle.pb_hp_model = inf(Np,1);

% Contact history fields used by update_particle.m.
Particle.delta_t = cell(Np, Np);
Particle.n_prev  = cell(Np, Np);
Particle.delta_t_wall = cell(Np, 4);

%% ========================================================================
% 7. Diagnostics storage
%% ========================================================================
max_steps_est = ceil(t_end / 1.0e-6) + 1000;

time_hist = nan(max_steps_est,1);
state_hist = nan(max_steps_est,1);
hp_hist = nan(max_steps_est,1);
timer_hist = nan(max_steps_est,1);
ta_hist = nan(max_steps_est,1);
Fcap_y_hist = nan(max_steps_est,1);
area_err_hist = nan(max_steps_est,1);
alpha_diff_hist = nan(max_steps_est,1);

%% ========================================================================
% 8. Visualization setup
%% ========================================================================
Xf = linspace(x_vec(1), x_vec(end), viz_up*N);
Yf = linspace(y_vec(1), y_vec(end), viz_up*M);
[XXf, YYf] = meshgrid(Xf, Yf);

fig = figure('Position', [80, 60, 900, 860], ...
             'Color', 'w', ...
             'Visible', 'on', ...
             'Name', 'Single Bubble - Single Particle Virtual Contact Test');

ax = axes('Parent', fig);
set(ax, 'YDir', 'normal');
axis(ax, 'equal');
axis(ax, [0 W 0 L]);
box(ax, 'on');
set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'Layer', 'top');
hold(ax, 'on');

colormap(ax, parula(256));
cb = colorbar(ax, 'Location', 'eastoutside');
cb.Label.String = 'Speed magnitude |u| (m/s)';

%% ========================================================================
% 9. Main loop
%% ========================================================================
timer = 0.0;
n = 0;

while true

    if timer >= t_end
        break;
    end

    n = n + 1;

    %% ---------------- timestep ----------------
    dt = compute_dt_cfl_liquid_only(Fluid, Particle, grid, XX, YY);
    dt = min(dt, t_end - timer);

    grid.dt = dt;
    timer = timer + dt;

    % Keep induction clips consistent with variable dt.
    grid.induction_t_min = 2 * grid.dt;
    grid.induction_t_max = 18 * grid.dt;

    %% ---------------- particle-bubble state ----------------
    if grid.enable_particle_bubble_state
        Particle = update_particle_bubble_state(Fluid, Particle, grid);
    end

    %% ---------------- real VOF advection ----------------
    % Do not use old update_alpha.
    alpha_real_new = update_advection(Fluid, grid);
    alpha_real_new = apply_boundary_conditions('vof', grid, alpha_real_new);

    %% ---------------- virtual extension alpha ----------------
    alpha_ext_new = alpha_real_new;

    if grid.enable_virtual_extension
        alpha_ext_new = update_extension(alpha_real_new, Particle, grid);
    end

    alpha_ext_new = apply_boundary_conditions('vof', grid, alpha_ext_new);

    % Store before prediction and Capillary.
    Fluid.alpha = alpha_real_new;
    Fluid.alpha_ext = alpha_ext_new;

    %% ---------------- capillary force on particle ----------------
    [F_cap_x, F_cap_y, T_cap] = Capillary(Fluid, Particle, grid);

    %% ---------------- fluid prediction ----------------
    [u_pred, v_pred] = prediction(Fluid, grid);

    %% ---------------- IBM / interaction ----------------
    [f_x, f_y] = Interaction(u_pred, v_pred, Particle, grid);

    %% ---------------- update particle ----------------
    Particle = update_particle(f_x, f_y, F_cap_x, F_cap_y, T_cap, ...
                               Particle, Fluid, grid, Wall);

    %% ---------------- update fluid / projection ----------------
    Fluid = update_fluid(u_pred, v_pred, ...
                         alpha_real_new, alpha_ext_new, ...
                         f_x, f_y, Fluid, grid);

    %% ---------------- diagnostics ----------------
    diag = compute_particle_bubble_diagnostics(Fluid, Particle, bub, grid, bubble_area0);

    time_hist(n) = timer;
    state_hist(n) = Particle.pb_state(1);
    hp_hist(n) = Particle.pb_hp_model(1);
    timer_hist(n) = Particle.induction_timer(1);
    ta_hist(n) = Particle.ta_induction(1);
    Fcap_y_hist(n) = F_cap_y(1);
    area_err_hist(n) = diag.area_err;
    alpha_diff_hist(n) = max(abs(Fluid.alpha_ext(:) - Fluid.alpha(:)));

    if mod(n, diag_every_nsteps) == 0 || n == 1
        fprintf(['step=%5d, t=%.4e, dt=%.2e, state=%d, ', ...
                 'hp_model=%.2e (%.2fh), timer=%.2e, ta=%.2e, ', ...
                 'area_err=%.2e, max|a_ext-a|=%.2e\n'], ...
                 n, timer, dt, Particle.pb_state(1), ...
                 Particle.pb_hp_model(1), Particle.pb_hp_model(1)/h, ...
                 Particle.induction_timer(1), Particle.ta_induction(1), ...
                 diag.area_err, alpha_diff_hist(n));
    end

    %% ---------------- visualization ----------------
    if mod(n, plot_every_nsteps) == 0 || n == 1
        cla(ax);
        hold(ax, 'on');
        set(ax, 'Color', 'w');

        u_cc = 0.5 * (Fluid.u(start:endy, start:endx) + ...
                      Fluid.u(start:endy, start+1:endx+1));

        v_cc = 0.5 * (Fluid.v(start:endy, start:endx) + ...
                      Fluid.v(start+1:endy+1, start:endx));

        speed = sqrt(u_cc.^2 + v_cc.^2);
        aplot = Fluid.alpha(start:endy, start:endx);

        speed_sm = smooth3x3(speed);
        aplot_sm = smooth3x3(aplot);

        speed_f = interp2(XX, YY, speed_sm, XXf, YYf, 'linear', 0.0);
        aplot_f = interp2(XX, YY, aplot_sm, XXf, YYf, 'linear', 1.0);

        u_f = interp2(XX, YY, u_cc, XXf, YYf, 'linear', 0.0);
        v_f = interp2(XX, YY, v_cc, XXf, YYf, 'linear', 0.0);

        % Speed magnitude contour.
        contourf(ax, XXf, YYf, speed_f, 24, 'LineStyle', 'none');
        colormap(ax, parula(256));
        caxis(ax, [0 speed_caxis_max]);

        % Sparse streamlines.
        sx = linspace(0.10*W, 0.90*W, 3);
        sy = linspace(0.10*L, 0.92*L, 4);
        [SX, SY] = meshgrid(sx, sy);

        try
            hs = streamslice(ax, XXf, YYf, u_f, v_f, 2);
            set(hs, 'Color', [0.05 0.05 0.05], 'LineWidth', 0.8);
        catch
            % Streamline can fail if velocity is almost exactly zero.
        end

        % Bubble interface alpha=0.5.
        try
            contour(ax, XXf, YYf, aplot_f, [0.5 0.5], ...
                    'k', 'LineWidth', 1.8);
        catch
        end

        % Particle real boundary.
        th = linspace(0, 2*pi, 240);

        xp = Particle.x_c(1) + Particle.r(1) * cos(th);
        yp = Particle.y_c(1) + Particle.r(1) * sin(th);

        patch(ax, xp, yp, [0.92 0.28 0.20], ...
              'EdgeColor', 'k', 'LineWidth', 1.2);

        % Virtual boundary.
        rv = Particle.r(1) + grid.virtual_film_thickness;
        xv = Particle.x_c(1) + rv * cos(th);
        yv = Particle.y_c(1) + rv * sin(th);

        plot(ax, xv, yv, ':', 'Color', [0.0 0.0 0.0], 'LineWidth', 1.2);

        % Initial bubble circle reference.
        xb0 = bub.x + bub.r * cos(th);
        yb0 = bub.y + bub.r * sin(th);
        plot(ax, xb0, yb0, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0);

        % Contact point from state detector.
        if isfield(Particle, 'pb_contact_x') && isfinite(Particle.pb_contact_x(1))
            plot(ax, Particle.pb_contact_x(1), Particle.pb_contact_y(1), ...
                 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
        end

        % Capillary force arrow.
        Fmag = sqrt(F_cap_x(1)^2 + F_cap_y(1)^2);
        if Fmag > 1e-14
            force_scale = 2.2 * Particle.r(1) / Fmag;
            quiver(ax, Particle.x_c(1), Particle.y_c(1), ...
                   force_scale * F_cap_x(1), force_scale * F_cap_y(1), ...
                   0, 'Color', [0.85 0.05 0.02], ...
                   'LineWidth', 2.0, 'MaxHeadSize', 1.8);
        end

        axis(ax, 'equal');
        axis(ax, [0 W 0 L]);
        box(ax, 'on');

        title(ax, 'Single Bubble - Single Particle Virtual Contact Test', ...
              'FontWeight', 'bold', 'FontSize', 13);

        xlabel(ax, 'x (m)');
        ylabel(ax, 'y (m)');

        state_text = state_name(Particle.pb_state(1));

        txt1 = sprintf('t = %.4e s   step = %d   dt = %.2e s', timer, n, dt);
        txt2 = sprintf('state = %d (%s)   active = %d', ...
                       Particle.pb_state(1), state_text, Particle.virtual_contact_active(1));
        txt3 = sprintf('h_p^{model} = %.2e m = %.2f h', ...
                       Particle.pb_hp_model(1), Particle.pb_hp_model(1)/h);
        txt4 = sprintf('timer = %.2e s   t_a = %.2e s', ...
                       Particle.induction_timer(1), Particle.ta_induction(1));
        txt5 = sprintf('Fcap_y = %.2e   area err = %.2e   max|a_e-a| = %.2e', ...
                       F_cap_y(1), diag.area_err, alpha_diff_hist(n));

        text(ax, 0.02, 0.965, sprintf('%s\n%s\n%s\n%s\n%s', ...
             txt1, txt2, txt3, txt4, txt5), ...
             'Units', 'normalized', ...
             'FontSize', 10.5, ...
             'FontWeight', 'bold', ...
             'Color', 'k', ...
             'BackgroundColor', 'w', ...
             'EdgeColor', 'k', ...
             'Margin', 3, ...
             'VerticalAlignment', 'top');

        drawnow;
    end
end

fprintf('Done. End time reached: %.6e s, steps = %d.\n', timer, n);

%% ========================================================================
% 10. Final diagnostics
%% ========================================================================
valid = isfinite(time_hist);
time_hist = time_hist(valid);
state_hist = state_hist(valid);
hp_hist = hp_hist(valid);
timer_hist = timer_hist(valid);
ta_hist = ta_hist(valid);
Fcap_y_hist = Fcap_y_hist(valid);
area_err_hist = area_err_hist(valid);
alpha_diff_hist = alpha_diff_hist(valid);

fig2 = figure('Position', [1020, 80, 680, 760], ...
              'Color', 'w', ...
              'Name', 'Diagnostics');

ax1 = subplot(5,1,1);
plot(ax1, time_hist, state_hist, 'LineWidth', 1.4);
ylabel(ax1, 'state');
ylim(ax1, [-0.2 3.2]);
box(ax1, 'on');

ax2 = subplot(5,1,2);
plot(ax2, time_hist, hp_hist/h, 'LineWidth', 1.4);
hold(ax2, 'on');
plot(ax2, [time_hist(1), time_hist(end)], [0, 0], 'k--');
ylabel(ax2, 'h_p/h');
box(ax2, 'on');

ax3 = subplot(5,1,3);
plot(ax3, time_hist, timer_hist, 'LineWidth', 1.4);
hold(ax3, 'on');
plot(ax3, time_hist, ta_hist, '--', 'LineWidth', 1.2);
ylabel(ax3, 'timer, t_a');
legend(ax3, 'timer', 't_a', 'Location', 'best');
box(ax3, 'on');

ax4 = subplot(5,1,4);
plot(ax4, time_hist, Fcap_y_hist, 'LineWidth', 1.4);
ylabel(ax4, 'F_{cap,y}');
box(ax4, 'on');

ax5 = subplot(5,1,5);
plot(ax5, time_hist, alpha_diff_hist, 'LineWidth', 1.4);
xlabel(ax5, 'time (s)');
ylabel(ax5, 'max|a_e-a|');
box(ax5, 'on');

%% ========================================================================
% Local helper functions
%% ========================================================================

function dt_val = compute_dt_cfl_liquid_only(Fluid, Particle, grid, XX, YY)

    start = grid.start;
    endx  = grid.endx;
    endy  = grid.endy;
    h     = grid.h;

    CFL_fluid = 0.35;
    CFL_part  = 0.25;
    Csig      = 0.20;
    Cvisc     = 0.25;

    dt_abs_max = 2.0e-5;
    dt_abs_min = 1.0e-9;
    max_growth = 1.25;

    u_cc = 0.5 * ...
        (Fluid.u(start:endy, start:endx) + ...
         Fluid.u(start:endy, start+1:endx+1));

    v_cc = 0.5 * ...
        (Fluid.v(start:endy, start:endx) + ...
         Fluid.v(start+1:endy+1, start:endx));

    speed = sqrt(u_cc.^2 + v_cc.^2);

    aplot = Fluid.alpha(start:endy, start:endx);
    mask_liq = (aplot > 1e-6);

    if isequal(size(XX), size(speed))
        XX_use = XX;
        YY_use = YY;
    else
        XX_use = XX(start:endy, start:endx);
        YY_use = YY(start:endy, start:endx);
    end

    mask_part = false(size(speed));

    if ~isempty(Particle) && isfield(Particle, 'r')
        for pp = 1:numel(Particle.r)
            mask_part = mask_part | ...
                ((XX_use - Particle.x_c(pp)).^2 + ...
                 (YY_use - Particle.y_c(pp)).^2 <= Particle.r(pp)^2);
        end
    end

    vals = speed(mask_liq & (~mask_part));

    if isempty(vals)
        Umax = 0.0;
    else
        Umax = max(vals(:));
    end

    if Umax < 1e-12
        dt_cfl = dt_abs_max;
    else
        dt_cfl = CFL_fluid * h / Umax;
    end

    rho_eff = Fluid.rhol;
    sigma = Fluid.sigma;
    dt_sigma = Csig * sqrt(rho_eff * h^3 / max(sigma, 1e-12));

    mu_max = max([Fluid.mul, Fluid.mug]);
    dt_visc = Cvisc * Fluid.rhol * h^2 / max(mu_max, 1e-14);

    dt_part = inf;

    if ~isempty(Particle)
        Up_max = 0.0;

        if isfield(Particle, 'u')
            Up_max = max(Up_max, max(abs(Particle.u(:))));
        end

        if isfield(Particle, 'v')
            Up_max = max(Up_max, max(abs(Particle.v(:))));
        end

        if Up_max > 1e-12
            dt_part = CFL_part * h / Up_max;
        end
    end

    dt_val = min([dt_cfl, dt_sigma, dt_visc, dt_part, dt_abs_max]);

    if isfield(grid, 'dt') && grid.dt > 0
        dt_val = min(dt_val, max_growth * grid.dt);
    end

    dt_val = max(dt_val, dt_abs_min);
end

function diag = compute_particle_bubble_diagnostics(Fluid, Particle, bub, grid, bubble_area0)

    start = grid.start;
    endx  = grid.endx;
    endy  = grid.endy;
    h     = grid.h;
    g     = grid.ghostnum;

    alpha = Fluid.alpha;

    gas_frac = 1.0 - alpha(start:endy, start:endx);
    gas_frac = min(max(gas_frac, 0.0), 1.0);

    bubble_area = sum(gas_frac(:)) * h^2;
    area_err = (bubble_area - bubble_area0) / max(bubble_area0, 1e-30);

    xp = Particle.x_c(1);
    yp = Particle.y_c(1);
    rp = Particle.r(1);

    gap_alpha = inf;
    x_int_min = NaN;
    y_int_min = NaN;

    for j = start+1:endy-1
        y = (j - g - 0.5) * h;

        for i = start+1:endx-1
            a = alpha(j,i);

            if a > 0.05 && a < 0.95
                x = (i - g - 0.5) * h;

                d_to_particle_surface = hypot(x - xp, y - yp) - rp;

                if d_to_particle_surface < gap_alpha
                    gap_alpha = d_to_particle_surface;
                    x_int_min = x;
                    y_int_min = y;
                end
            end
        end
    end

    d_center = hypot(xp - bub.x, yp - bub.y);
    gap_geom = d_center - rp - bub.r;

    diag.gap_alpha = gap_alpha;
    diag.gap_geom = gap_geom;
    diag.x_int_min = x_int_min;
    diag.y_int_min = y_int_min;
    diag.bubble_area = bubble_area;
    diag.area_err = area_err;
end

function Aout = smooth3x3(Ain)

    K = [1 2 1; 2 4 2; 1 2 1] / 16;
    Aout = Ain;

    [mm, nn] = size(Ain);

    for jj = 2:mm-1
        for ii = 2:nn-1
            patch = Ain(jj-1:jj+1, ii-1:ii+1);
            Aout(jj,ii) = sum(sum(patch .* K));
        end
    end
end

function alpha_out = embed_bubbles_union(alpha_in, bub, grid, subN)
% Subcell-sampled circular gas bubble.
% alpha = 1 liquid, alpha = 0 gas.

    alpha_out = alpha_in;

    ghostnum = grid.ghostnum;
    start = grid.start;
    endx = grid.endx;
    endy = grid.endy;
    h = grid.h;

    dx = h / subN;
    dy = h / subN;

    nb = numel(bub.r);
    maxR = max(bub.r);

    for j = start:endy
        for i = start:endx

            xmin = (i - ghostnum - 1) * h;
            ymin = (j - ghostnum - 1) * h;

            cx = xmin + 0.5*h;
            cy = ymin + 0.5*h;

            if all(hypot(cx - bub.x(:), cy - bub.y(:)) > (maxR + 2*h))
                alpha_out(j,i) = 1.0;
                continue;
            end

            insideCount = 0;

            for py = 1:subN
                yy = ymin + (py - 0.5) * dy;

                for px = 1:subN
                    xx = xmin + (px - 0.5) * dx;

                    insideAny = false;

                    for kk = 1:nb
                        if ((xx - bub.x(kk))^2 + (yy - bub.y(kk))^2) <= bub.r(kk)^2
                            insideAny = true;
                            break;
                        end
                    end

                    if insideAny
                        insideCount = insideCount + 1;
                    end
                end
            end

            gasFrac = insideCount / (subN * subN);
            alpha_out(j,i) = 1.0 - gasFrac;
        end
    end
end

function alpha_new = fill_vof_ghost_neumann(alpha_in, grid)

    alpha_new = alpha_in;

    g = grid.ghostnum;
    start = grid.start;
    endx = grid.endx;
    endy = grid.endy;
    Ny = grid.Ny;
    Nx = grid.Nx;

    for i = 1:g
        alpha_new(:, i) = alpha_new(:, start);
    end

    for i = endx+1:Nx
        alpha_new(:, i) = alpha_new(:, endx);
    end

    for j = 1:g
        alpha_new(j, :) = alpha_new(start, :);
    end

    for j = endy+1:Ny
        alpha_new(j, :) = alpha_new(endy, :);
    end
end

function name = state_name(s)

    switch s
        case 0
            name = 'free';
        case 1
            name = 'collision';
        case 2
            name = 'induction';
        case 3
            name = 'attachment';
        otherwise
            name = 'unknown';
    end
end