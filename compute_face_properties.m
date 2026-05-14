function prop = compute_face_properties(Fluid, grid)
%% ========================================================================
% compute_face_properties
%
% Unified material properties for prediction and Poisson correction.
%
% alpha convention:
%   alpha = 1 liquid
%   alpha = 0 gas
%
% Returns:
%   prop.rho_cc
%   prop.mu_cc
%   prop.invrho_x   face-centered 1/rho for u-faces
%   prop.invrho_y   face-centered 1/rho for v-faces
%   prop.mu_x       face-centered mu for u-faces
%   prop.mu_y       face-centered mu for v-faces
%
% NOTE:
%   This version restores the arithmetic average of cell-centered 1/rho.
%   It is safer for the current balanced-force consistency of your code.
%% ========================================================================

start = grid.start;
endy  = grid.endy;
endx  = grid.endx;
Ny    = grid.Ny;
Nx    = grid.Nx;

eps_div = 1e-12;

alpha = Fluid.alpha;
alpha = apply_boundary_conditions('vof', grid, alpha);

rhol = Fluid.rhol;
rhog = Fluid.rhog;
mul  = Fluid.mul;
mug  = Fluid.mug;

rho_cc = rhol .* alpha + rhog .* (1 - alpha);
mu_cc  = mul  .* alpha + mug  .* (1 - alpha);

inv_rho_cc = 1 ./ (rho_cc + eps_div);

invrho_x = zeros(Ny+1, Nx+1);
invrho_y = zeros(Ny+1, Nx+1);
mu_x     = zeros(Ny+1, Nx+1);
mu_y     = zeros(Ny+1, Nx+1);

for j = start:endy+1
    for i = start:endx+1

        % u-face: between cell (j,i-1) and cell (j,i)
        invrho_x(j,i) = 0.5 * (inv_rho_cc(j,i-1) + inv_rho_cc(j,i));
        mu_x(j,i)     = 0.5 * (mu_cc(j,i-1)      + mu_cc(j,i));

        % v-face: between cell (j-1,i) and cell (j,i)
        invrho_y(j,i) = 0.5 * (inv_rho_cc(j-1,i) + inv_rho_cc(j,i));
        mu_y(j,i)     = 0.5 * (mu_cc(j-1,i)      + mu_cc(j,i));
    end
end

prop.rho_cc    = rho_cc;
prop.mu_cc     = mu_cc;
prop.invrho_x  = invrho_x;
prop.invrho_y  = invrho_y;
prop.mu_x      = mu_x;
prop.mu_y      = mu_y;

end