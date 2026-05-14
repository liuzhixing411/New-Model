function Fluid = update_fluid(u_pred, v_pred, alpha_new, varargin)
% =========================================================================
% update_fluid
%
% Incremental projection version with separated real alpha and alpha_ext.
%
% Main purpose:
%   Fluid.alpha     : real VOF field
%                     used for mass conservation, density, viscosity,
%                     Poisson variable coefficient, etc.
%
%   Fluid.alpha_ext : auxiliary capillary-geometry field
%                     used only by CSF / Capillary if corresponding
%                     switches are enabled.
%
% New recommended calling form:
%
%   Fluid = update_fluid(u_pred, v_pred, alpha_new, alpha_ext_new, ...
%                        f_x, f_y, Fluid, grid);
%
% Backward compatible calling form:
%
%   Fluid = update_fluid(u_pred, v_pred, alpha_new, ...
%                        f_x, f_y, Fluid, grid);
%
% Projection convention:
%   prediction(...) has already used old pressure p^n.
%   Poisson(...) solves pressure correction phi.
%
% Steps:
%   1) u_star = u_pred + f*dt
%   2) update Fluid.alpha and Fluid.alpha_ext before projection
%   3) solve pressure correction and project velocity
%   4) write u, v, alpha, alpha_ext, p back to Fluid
%
% IMPORTANT:
%   Poisson and compute_face_properties should use Fluid.alpha, not
%   Fluid.alpha_ext.
% =========================================================================

%% ------------------------------------------------------------------------
% Parse inputs
%% ------------------------------------------------------------------------
% New form:
%   varargin = {alpha_ext_new, f_x, f_y, Fluid, grid}
%
% Old form:
%   varargin = {f_x, f_y, Fluid, grid}
%% ------------------------------------------------------------------------

if numel(varargin) == 5

    alpha_ext_new = varargin{1};
    f_x           = varargin{2};
    f_y           = varargin{3};
    Fluid         = varargin{4};
    grid          = varargin{5};

elseif numel(varargin) == 4

    alpha_ext_new = [];
    f_x           = varargin{1};
    f_y           = varargin{2};
    Fluid         = varargin{3};
    grid          = varargin{4};

else
    error(['update_fluid: invalid input number. ', ...
           'Use either update_fluid(u_pred,v_pred,alpha_new,f_x,f_y,Fluid,grid) ', ...
           'or update_fluid(u_pred,v_pred,alpha_new,alpha_ext_new,f_x,f_y,Fluid,grid).']);
end

dt = grid.dt;

%% ------------------------------------------------------------------------
% Step 1: add IBM / interaction force to predicted velocity
%% ------------------------------------------------------------------------
u_star = u_pred + f_x * dt;
v_star = v_pred + f_y * dt;

%% ------------------------------------------------------------------------
% Step 2: update real alpha before projection
%
% Fluid.alpha is the real conservative VOF field.
% It must be updated before Poisson so that density / invrho are based on
% the new real alpha.
%% ------------------------------------------------------------------------
Fluid.alpha = alpha_new;
Fluid.alpha = apply_boundary_conditions('vof', grid, Fluid.alpha);

%% ------------------------------------------------------------------------
% Step 2b: update alpha_ext
%
% alpha_ext is only a capillary geometry field.
% If alpha_ext_new is not provided, fall back safely.
%% ------------------------------------------------------------------------
if isempty(alpha_ext_new)

    if isfield(Fluid, 'alpha_ext') && ~isempty(Fluid.alpha_ext)
        % Keep existing alpha_ext if it already exists.
        alpha_ext_new = Fluid.alpha_ext;
    else
        % Otherwise use real alpha.
        alpha_ext_new = Fluid.alpha;
    end
end

Fluid.alpha_ext = alpha_ext_new;
Fluid.alpha_ext = apply_boundary_conditions('vof', grid, Fluid.alpha_ext);

%% ------------------------------------------------------------------------
% Step 3: projection
%
% Poisson should use Fluid.alpha internally through compute_face_properties.
% It should NOT use Fluid.alpha_ext.
%% ------------------------------------------------------------------------
[u_new, v_new, p_new] = Poisson(Fluid, grid, u_star, v_star);

%% ------------------------------------------------------------------------
% Step 4: write back fields
%% ------------------------------------------------------------------------
Fluid.u = u_new;
Fluid.v = v_new;
Fluid.p = p_new;

% Re-apply boundary conditions.
Fluid.alpha     = apply_boundary_conditions('vof', grid, Fluid.alpha);
Fluid.alpha_ext = apply_boundary_conditions('vof', grid, Fluid.alpha_ext);

[Fluid.u, Fluid.v] = apply_boundary_conditions('velocity', grid, Fluid.u, Fluid.v);
Fluid.p = apply_boundary_conditions('pressure', grid, Fluid.p);

end