function [f_x,f_y]=Interaction(u_pred,v_pred,Particle,grid)

Ny=grid.Ny;
Nx=grid.Nx;
h=grid.h;
dt=grid.dt;
ghostnum=grid.ghostnum;
% Initial Force
f_x = zeros(Ny+1, Nx+1);
f_y = zeros(Ny+1, Nx+1);

for p=1:length(Particle.r)
    % Particle Parameter
    x_c = Particle.x_c(p);
    y_c = Particle.y_c(p);
    r = Particle.r(p);
    U_p = Particle.u(p);
    V_p = Particle.v(p);
    omega = Particle.omega(p);


    % Dirac delta function
    delta_kernel = @(r_val) ...
        (abs(r_val) >= 0 & abs(r_val) <= 1) .* (1/8) .* (3 - 2*abs(r_val) + sqrt(1 + 4*abs(r_val) - 4*abs(r_val).^2)) + ...
        (abs(r_val) > 1 & abs(r_val) <= 2) .* (1/8) .* (5 - 2*abs(r_val) - sqrt(-7 + 12*abs(r_val) - 4*abs(r_val).^2)) + ...
        (abs(r_val) > 2) .* 0;

    dirac = @(dx,dy) (1/h^2) * delta_kernel(dx/h) * delta_kernel(dy/h);

    % Generate Lagrangian Points
    N_lag = 100;  % Number of Lagrangian Points
    theta = linspace(0, 2*pi, N_lag+1);
    theta = theta(1:end-1);  % Remove repetitive points

    % Lagrangian Points Coordinates
    x_lag = x_c + r * cos(theta);
    y_lag = y_c + r * sin(theta);

    % Lagrangian Arc Length
    dS = 2*pi*r / N_lag*h;



    %% Project u_pred and v_pred to Lagrangian Points
    u_hat_lag = zeros(1, N_lag);
    v_hat_lag = zeros(1, N_lag);

    for k = 1:N_lag
        x_k = x_lag(k);
        y_k = y_lag(k);

        % Relevant Eulerian Points (From -2 to +2)
        for jp = -2:2
            for ip = -2:2
                % u-face
                j_u = round((y_k)/h+0.5) + ghostnum + jp;
                i_u = round((x_k)/h+1) + ghostnum + ip;

                if j_u >= 1 && j_u <= Ny+1 && i_u >= 1 && i_u <= Nx+1
                    x_u = (i_u - ghostnum - 1) * h;
                    y_u = (j_u - ghostnum - 0.5) * h;
                    dx = x_k - x_u;
                    dy = y_k - y_u;
                    delta_val = dirac(dx, dy);
                    u_hat_lag(k) = u_hat_lag(k) + u_pred(j_u,i_u) * delta_val * h^2;
                end

                % v-face
                j_v = round((y_k)/h+1) + ghostnum + jp;
                i_v = round((x_k)/h+0.5) + ghostnum + ip;

                if j_v >= 1 && j_v <= Ny+1 && i_v >= 1 && i_v <= Nx+1
                    x_v = (i_v - ghostnum-0.5) * h;
                    y_v = (j_v - ghostnum-1) * h;
                    dx = x_k - x_v;
                    dy = y_k - y_v;
                    delta_val = dirac(dx, dy);
                    v_hat_lag(k) = v_hat_lag(k) + v_pred(j_v,i_v) * delta_val * h^2;
                end
            end
        end
    end

    %% Step2：Calculate expected velocity on Lagrangian
    u_L_lag = zeros(1, N_lag);
    v_L_lag = zeros(1, N_lag);

    for k = 1:N_lag
        x_k = x_lag(k);
        y_k = y_lag(k);

        % 期望速度：u_L = U - ω*(y - y_c), v_L = V + ω*(x - x_c)
        u_L_lag(k) = U_p - omega * (y_k - y_c);
        v_L_lag(k) = V_p + omega * (x_k - x_c);
    end

    %% 步骤3：计算Lagrangian点上的力
    F_x_lag = (u_L_lag - u_hat_lag) / dt;
    F_y_lag = (v_L_lag - v_hat_lag) / dt;

    %% 步骤4：将Lagrangian点上的力分布到Eulerian网格
    for k = 1:N_lag
        x_k = x_lag(k);
        y_k = y_lag(k);

        % 遍历周围的Eulerian点
        for jp = -2:2
            for ip = -2:2
                % 分布f_x到u-faces
                j_u = round((y_k)/h+0.5 ) + ghostnum + jp;
                i_u = round((x_k)/h+1 ) + ghostnum + ip;

                if j_u >= 1 && j_u <= Ny+1 && i_u >= 1 && i_u <= Nx+1
                    x_u = (i_u - ghostnum - 1) * h;
                    y_u = (j_u - ghostnum - 0.5) * h;
                    dx = x_k - x_u;
                    dy = y_k - y_u;
                    delta_val = dirac(dx, dy);
                    f_x(j_u,i_u) = f_x(j_u,i_u) + F_x_lag(k) * delta_val * dS;
                end

                % 分布f_y到v-faces
                j_v = round((y_k)/h+1) + ghostnum + jp;
                i_v = round((x_k)/h+0.5) + ghostnum + ip;

                if j_v >= 1 && j_v <= Ny+1 && i_v >= 1 && i_v <= Nx+1
                    x_v = (i_v - ghostnum-0.5) * h;
                    y_v = (j_v - ghostnum-1) * h;
                    dx = x_k - x_v;
                    dy = y_k - y_v;
                    delta_val = dirac(dx, dy);
                    f_y(j_v,i_v) = f_y(j_v,i_v) + F_y_lag(k) * delta_val * dS;
                end
            end
        end
    end
end
end