function Particle = update_particle(f_x, f_y, F_cap_x, F_cap_y, T_cap, Particle, Fluid, grid, Wall)
% UPDATE_PARTICLE
% 更新粒子运动，包含：
% 1. 流体动力 (IBM 宏观力)
% 2. 润滑力修正
% 3. 软球碰撞模型
% 4. 子步时间积分
%
% 本版本修正：
%   A) F_hydro 的索引映射加上 ghost 偏移
%   B) F_hydro 只在粒子附近圆形区域积分
%   C) 浮力改为局部有效密度 rho_eff(alpha)
%
% 参考: Costa et al. 2015

%% 1. 物理参数
N_col   = 8;
e_n_d   = 0.97;
e_t_d   = 0.39;
mu_c    = 0.15;
K_gyr   = 1/2;

% 润滑力参数
eps_dx    = 0.05;
eps_sigma = 0.001;
mu_fluid  = Fluid.mul;

% 子步积分
N_sub  = 10;
dt_sub = grid.dt / N_sub;

%% 2. 准备工作
Np = length(Particle.r);
if ~isfield(Particle, 'delta_t') || isempty(Particle.delta_t)
    Particle.delta_t = cell(Np, Np);
    Particle.n_prev  = cell(Np, Np);
    Particle.delta_t_wall = cell(Np, 4);
end

dV = grid.h * grid.h;
rho_mix = Fluid.rhol .* Fluid.alpha + Fluid.rhog .* (1 - Fluid.alpha);
g_vec = [0; -9.81];

%% 3. 计算流体宏观力 F_hydro
F_hydro_x = zeros(Np, 1);
F_hydro_y = zeros(Np, 1);
T_hydro   = zeros(Np, 1);

for p = 1:Np
    x_c = Particle.x_c(p);
    y_c = Particle.y_c(p);
    r_p = Particle.r(p);

    % --- 正确的物理坐标 -> 数组索引映射 ---
    ic = round(x_c / grid.h + grid.ghostnum + 0.5);
    jc = round(y_c / grid.h + grid.ghostnum + 0.5);

    Rbox = ceil(r_p / grid.h) + 3;

    start_i = max(grid.start, ic - Rbox);
    end_i   = min(grid.endx+1, ic + Rbox);
    start_j = max(grid.start, jc - Rbox);
    end_j   = min(grid.endy+1, jc + Rbox);

    fx_sum = 0.0;
    fy_sum = 0.0;
    t_sum  = 0.0;

    for j = start_j:end_j
        for i = start_i:end_i

            % ---- U-face contribution ----
            if i > 1 && i <= size(f_x,2) && j >= 1 && j <= size(f_x,1)
                % u-face physical coordinates
                xg = (i - grid.ghostnum - 1.0) * grid.h;
                yg = (j - grid.ghostnum - 0.5) * grid.h;

                % only integrate near particle
                if (xg - x_c)^2 + (yg - y_c)^2 <= (1.2*r_p)^2
                    rho_u = 0.5 * (rho_mix(min(j, size(rho_mix,1)), min(i, size(rho_mix,2))) + ...
                        rho_mix(min(j, size(rho_mix,1)), max(i-1, 1)));

                    force_val = -(rho_u * f_x(j,i) * dV);
                    fx_sum = fx_sum + force_val;

                    % torque arm about particle center
                    dist_y = yg - y_c;
                    t_sum = t_sum - dist_y * force_val;
                end
            end

            % ---- V-face contribution ----
            if j > 1 && j <= size(f_y,1) && i >= 1 && i <= size(f_y,2)
                % v-face physical coordinates
                xg = (i - grid.ghostnum - 0.5) * grid.h;
                yg = (j - grid.ghostnum - 1.0) * grid.h;

                if (xg - x_c)^2 + (yg - y_c)^2 <= (1.2*r_p)^2
                    rho_v = 0.5 * (rho_mix(min(j, size(rho_mix,1)), min(i, size(rho_mix,2))) + ...
                        rho_mix(max(j-1, 1), min(i, size(rho_mix,2))));

                    force_val = -(rho_v * f_y(j,i) * dV);
                    fy_sum = fy_sum + force_val;

                    dist_x = xg - x_c;
                    t_sum = t_sum + dist_x * force_val;
                end
            end

        end
    end

    F_hydro_x(p) = fx_sum;
    F_hydro_y(p) = fy_sum;
    T_hydro(p)   = t_sum;
end

%% 4. 子步积分
for step = 1:N_sub
    for p = 1:Np
        m = Particle.m(p);
        I = Particle.I(p);
        r_p = Particle.r(p);
        pos_p = [Particle.x_c(p); Particle.y_c(p)];
        vel_p = [Particle.u(p); Particle.v(p)];
        omega_p = Particle.omega(p);

        F_cont_x = 0.0;
        F_cont_y = 0.0;
        T_cont   = 0.0;

        % -------------------------------------------------------
        % A. 粒子-粒子 交互
        % -------------------------------------------------------
        for q = 1:Np
            if p == q
                continue;
            end

            pos_q = [Particle.x_c(q); Particle.y_c(q)];
            vel_q = [Particle.u(q); Particle.v(q)];
            omega_q = Particle.omega(q);
            r_q = Particle.r(q);
            m_q = Particle.m(q);

            rel_pos = pos_q - pos_p;
            dist = norm(rel_pos);
            if dist == 0
                continue;
            end

            n_vec = rel_pos / dist;   % p -> q
            gap = dist - (r_p + r_q);
            eps = gap / max(r_p, 1e-12);

            % contact-point relative velocity
            v_surf_p = vel_p + cross_product_2d(omega_p, r_p * n_vec);
            v_surf_q = vel_q + cross_product_2d(omega_q, -r_q * n_vec);
            u_ij = v_surf_p - v_surf_q;

            u_n_val = dot(u_ij, n_vec);
            u_n_vec = u_n_val * n_vec;
            u_t_vec = u_ij - u_n_vec;

            % ---- 1) lubrication ----
            % approaching if u_n_val > 0 under your convention
            if gap > 0 && gap < eps_dx * grid.h && u_n_val > 0
                lambda_eps = lubrication_factor_pp(eps, eps_sigma, eps_dx * grid.h / max(r_p,1e-12));
                F_lub_val = -6 * pi * mu_fluid * r_p * u_n_val * lambda_eps;
                F_lub = F_lub_val * n_vec;

                F_cont_x = F_cont_x + F_lub(1);
                F_cont_y = F_cont_y + F_lub(2);
            end

            % ---- 2) collision ----
            overlap = -gap;
            if overlap > 0
                m_eff = (m * m_q) / (m + m_q);
                [kn, eta_n] = calc_coeffs(m_eff, e_n_d, N_col * grid.dt);

                F_n = -kn * overlap * n_vec - eta_n * u_n_vec;

                m_eff_t = m_eff / (1 + 1/K_gyr);
                [kt, eta_t] = calc_coeffs(m_eff_t, e_t_d, N_col * grid.dt);

                idxs = sort([p, q]);
                idx1 = idxs(1);
                idx2 = idxs(2);

                if isempty(Particle.delta_t{idx1, idx2})
                    Particle.delta_t{idx1, idx2} = [0;0];
                    Particle.n_prev{idx1, idx2} = n_vec;
                end

                n_prev = Particle.n_prev{idx1, idx2};
                delta_t_old = Particle.delta_t{idx1, idx2};

                t_vec_curr = [-n_vec(2); n_vec(1)];
                t_vec_prev = [-n_prev(2); n_prev(1)];
                delta_scalar = dot(delta_t_old, t_vec_prev);
                delta_t_rot = delta_scalar * t_vec_curr;

                delta_t_star = delta_t_rot + u_t_vec * dt_sub;

                F_t_trial = -kt * delta_t_star - eta_t * u_t_vec;
                Fn_mag = norm(F_n);
                Ft_mag = norm(F_t_trial);

                if Ft_mag <= mu_c * Fn_mag
                    F_t = F_t_trial;
                    Particle.delta_t{idx1, idx2} = delta_t_star;
                else
                    t_dir = F_t_trial / (Ft_mag + 1e-12);
                    F_t = mu_c * Fn_mag * t_dir;
                    Particle.delta_t{idx1, idx2} = (-F_t - eta_t * u_t_vec) / kt;
                end
                Particle.n_prev{idx1, idx2} = n_vec;

                F_total_ij = F_n + F_t;
                F_cont_x = F_cont_x + F_total_ij(1);
                F_cont_y = F_cont_y + F_total_ij(2);

                T_cont = T_cont + (r_p * n_vec(1) * F_total_ij(2) - r_p * n_vec(2) * F_total_ij(1));
            else
                idxs = sort([p, q]);
                idx1 = idxs(1);
                idx2 = idxs(2);

                if ~isempty(Particle.delta_t) && numel(Particle.delta_t) >= (idx1-1)*Np + idx2
                    Particle.delta_t{idx1, idx2} = [];
                end
                if ~isempty(Particle.n_prev) && numel(Particle.n_prev) >= (idx1-1)*Np + idx2
                    Particle.n_prev{idx1, idx2} = [];
                end
            end
        end

        % -------------------------------------------------------
        % B. 粒子-墙 交互
        % -------------------------------------------------------
        normals = {[1;0], [-1;0], [0;1], [0;-1]};
        side_names = {'left', 'right', 'down', 'up'};

        for w = 1:4

            side_name = side_names{w};
            side_bc = get_wall_bc_type(grid, side_name);

            % 只有 non-slip / no-slip 才视为真实固壁
            % free-slip 和 periodic 均不计算墙面润滑和墙面碰撞
            if ~strcmp(side_bc, 'non-slip')
                if isfield(Particle, 'delta_t_wall') && numel(Particle.delta_t_wall) >= (p-1)*4 + w
                    Particle.delta_t_wall{p, w} = [];
                end
                continue;
            end

            n_w = normals{w};

            switch w
                case 1 % left
                    dist_in = pos_p(1) - Wall.left;
                case 2 % right
                    dist_in = Wall.right - pos_p(1);
                case 3 % bottom
                    dist_in = pos_p(2) - Wall.down;
                case 4 % top
                    dist_in = Wall.up - pos_p(2);
            end

            gap = dist_in - r_p;
            eps = gap / max(r_p, 1e-12);

            v_surf_p = vel_p + cross_product_2d(omega_p, r_p * (-n_w));
            u_pw = v_surf_p;

            u_n = max(0, -dot(u_pw, n_w));

            % ---- lubrication: only for solid wall ----
            if gap > 0 && gap < eps_dx * grid.h && u_n > 0
                lambda_eps = lubrication_factor_pw( ...
                    max(eps, eps_sigma), ...
                    eps_sigma, ...
                    eps_dx * grid.h / max(r_p,1e-12));

                F_lub_mag = 6 * pi * mu_fluid * r_p * u_n * lambda_eps;
                F_lub = F_lub_mag * n_w;

                F_cont_x = F_cont_x + F_lub(1);
                F_cont_y = F_cont_y + F_lub(2);
            end

            % ---- contact: only for solid wall ----
            overlap = -gap;
            if overlap > 0
                m_eff = m;
                [kn, eta_n] = calc_coeffs(m_eff, e_n_d, N_col * grid.dt);

                u_n_vec = (dot(u_pw, n_w)) * n_w;
                u_t_vec = u_pw - u_n_vec;

                F_n = kn * overlap * n_w - eta_n * u_n_vec;

                m_eff_t = m_eff / (1 + 1/K_gyr);
                [kt, eta_t] = calc_coeffs(m_eff_t, e_t_d, N_col * grid.dt);

                if isempty(Particle.delta_t_wall{p, w})
                    Particle.delta_t_wall{p, w} = [0;0];
                end

                delta_t_star = Particle.delta_t_wall{p, w} + u_t_vec * dt_sub;

                F_t_trial = -kt * delta_t_star - eta_t * u_t_vec;

                if norm(F_t_trial) <= mu_c * norm(F_n)
                    F_t = F_t_trial;
                    Particle.delta_t_wall{p, w} = delta_t_star;
                else
                    t_dir = F_t_trial / (norm(F_t_trial)+1e-12);
                    F_t = mu_c * norm(F_n) * t_dir;
                    Particle.delta_t_wall{p, w} = (-F_t - eta_t * u_t_vec) / kt;
                end

                F_total = F_n + F_t;

                F_cont_x = F_cont_x + F_total(1);
                F_cont_y = F_cont_y + F_total(2);

                vec_arm = -r_p * n_w;
                T_cont = T_cont + ...
                    (vec_arm(1)*F_total(2) - vec_arm(2)*F_total(1));

            else
                Particle.delta_t_wall{p, w} = [];
            end
        end
        % -------------------------------------------------------
        % C. 体力：局部有效密度浮力修正
        % -------------------------------------------------------
        ic_loc = min(max(round(Particle.x_c(p)/grid.h + grid.ghostnum + 0.5), grid.start), grid.endx);
        jc_loc = min(max(round(Particle.y_c(p)/grid.h + grid.ghostnum + 0.5), grid.start), grid.endy);

        alpha_c = Fluid.alpha(jc_loc, ic_loc);
        rho_f_approx = alpha_c * Fluid.rhol + (1 - alpha_c) * Fluid.rhog;

        F_buoy_x = 0.0;
        F_buoy_y = (Particle.rho(p) - rho_f_approx) * Particle.V(p) * g_vec(2);

        % -------------------------------------------------------
        % D. 合力与积分
        % -------------------------------------------------------
        Total_Fx = F_hydro_x(p) + F_buoy_x + F_cont_x + F_cap_x(p);
        Total_Fy = F_hydro_y(p) + F_buoy_y + F_cont_y + F_cap_y(p);
        Total_T  = T_hydro(p)   + T_cont   + T_cap(p);

        % 调试时可打开
        % if p == 1 && step == 1
        %     fprintf('Fhydro_y=%.4e  Fbuoy_y=%.4e  Fcont_y=%.4e  Fcap_y=%.4e\n', ...
        %         F_hydro_y(p), F_buoy_y, F_cont_y, F_cap_y(p));
        % end

        Particle.u(p) = Particle.u(p) + (Total_Fx / m) * dt_sub;
        Particle.v(p) = Particle.v(p) + (Total_Fy / m) * dt_sub;
        Particle.omega(p) = Particle.omega(p) + (Total_T / I) * dt_sub;

        Particle.x_c(p) = Particle.x_c(p) + Particle.u(p) * dt_sub;
        Particle.y_c(p) = Particle.y_c(p) + Particle.v(p) * dt_sub;
    end
end

end

%% ============================================================
% 辅助函数
%% ============================================================

function [kn, eta] = calc_coeffs(m_eff, e, T_col)
if T_col < 1e-9
    T_col = 1e-9;
end
ln_e = log(e);
factor = sqrt(pi^2 + ln_e^2);
kn = m_eff * (factor / T_col)^2;
eta = -2 * m_eff * ln_e / T_col;
end

function res = cross_product_2d(omega_scalar, r_vec)
res = [-omega_scalar * r_vec(2); omega_scalar * r_vec(1)];
end

function lambda = lubrication_factor_pp(eps, eps_sigma, eps_dx)
val_eps = max(eps, eps_sigma);
lambda_raw = 1/(2*val_eps);
lambda_cutoff = 1/(2*eps_dx);

lambda = lambda_raw - lambda_cutoff;
if lambda < 0
    lambda = 0;
end
end

function lambda = lubrication_factor_pw(eps, eps_sigma, eps_dx)
val_eps = max(eps, eps_sigma);
lambda_raw = 1/val_eps;
lambda_cutoff = 1/eps_dx;

lambda = lambda_raw - lambda_cutoff;
if lambda < 0
    lambda = 0;
end
end

function side_bc = get_wall_bc_type(grid, side_name)

if ~isfield(grid, 'bc')
    side_bc = 'non-slip';
    return;
end

if ~isfield(grid.bc, side_name)
    side_bc = 'non-slip';
    return;
end

side_bc = lower(strtrim(grid.bc.(side_name)));
side_bc = strrep(side_bc, '_', '-');

switch side_bc
    case {'no-slip','noslip','non-slip','nonslip','solid'}
        side_bc = 'non-slip';

    case {'free-slip','freeslip','slip'}
        side_bc = 'free-slip';

    case {'periodic','perodic'}
        side_bc = 'periodic';

    otherwise
        error('Unknown wall boundary type: %s', side_bc);
end

end