function varargout = apply_boundary_conditions(type, grid, varargin)
% APPLY_BOUNDARY_CONDITIONS
%
% Unified boundary-condition wrapper.
%
% Usage:
%   alpha = apply_boundary_conditions('scalar', grid, alpha);
%   p     = apply_boundary_conditions('pressure', grid, p);
%   [u,v] = apply_boundary_conditions('velocity', grid, u, v);
%
% Supported boundary types in grid.bc:
%   'non-slip' / 'no-slip'
%   'free-slip'
%   'periodic'

switch lower(type)
    case {'scalar','pressure','vof'}
        A = varargin{1};
        A = apply_scalar_bc(A, grid);
        varargout{1} = A;

    case 'velocity'
        u = varargin{1};
        v = varargin{2};
        [u, v] = apply_velocity_bc(u, v, grid);
        varargout{1} = u;
        varargout{2} = v;

    otherwise
        error('Unknown boundary condition field type: %s', type);
end

end

%% ========================================================================
% Scalar boundary condition: alpha, pressure, density-like fields
%% ========================================================================
function A = apply_scalar_bc(A, grid)

g     = grid.ghostnum;
start = grid.start;
endx  = grid.endx;
endy  = grid.endy;
Ny    = size(A,1);
Nx    = size(A,2);

bc = get_bc(grid);

% ---------------- left / right ----------------
if is_periodic_pair(bc.left, bc.right)

    for k = 1:g
        A(:, start-k) = A(:, endx-k+1);
        A(:, endx+k)  = A(:, start+k-1);
    end

else
    % left
    switch normalize_bc(bc.left)
        case {'non-slip','free-slip'}
            for k = 1:g
                A(:, start-k) = A(:, start);
            end
        otherwise
            error('Unsupported left scalar BC: %s', bc.left);
    end

    % right
    switch normalize_bc(bc.right)
        case {'non-slip','free-slip'}
            for k = 1:g
                A(:, endx+k) = A(:, endx);
            end
        otherwise
            error('Unsupported right scalar BC: %s', bc.right);
    end
end

% ---------------- down / up ----------------
if is_periodic_pair(bc.down, bc.up)

    for k = 1:g
        A(start-k, :) = A(endy-k+1, :);
        A(endy+k, :)  = A(start+k-1, :);
    end

else
    % down
    switch normalize_bc(bc.down)
        case {'non-slip','free-slip'}
            for k = 1:g
                A(start-k, :) = A(start, :);
            end
        otherwise
            error('Unsupported down scalar BC: %s', bc.down);
    end

    % up
    switch normalize_bc(bc.up)
        case {'non-slip','free-slip'}
            for k = 1:g
                A(endy+k, :) = A(endy, :);
            end
        otherwise
            error('Unsupported up scalar BC: %s', bc.up);
    end
end

% 防止角点由于先后覆盖出现 NaN
A = fill_scalar_corners(A, start, endx, endy, Ny, Nx);

end

%% ========================================================================
% Velocity boundary condition on staggered grid
%% ========================================================================
function [u, v] = apply_velocity_bc(u, v, grid)

g     = grid.ghostnum;
start = grid.start;
endx  = grid.endx;
endy  = grid.endy;

Ny_u = size(u,1);
Nx_u = size(u,2);
Ny_v = size(v,1);
Nx_v = size(v,2);

bc = get_bc(grid);

%% ---------------- x-direction boundaries: left/right ----------------
if is_periodic_pair(bc.left, bc.right)

    % u-face periodic wrap
    for k = 1:g
        u(:, start-k)   = u(:, endx+1-k);
        u(:, endx+1+k) = u(:, start+k);
    end

    % v-face periodic wrap
    for k = 1:g
        v(:, start-k) = v(:, endx-k+1);
        v(:, endx+k)  = v(:, start+k-1);
    end

else
    % ---------------- left boundary ----------------
    left_bc = normalize_bc(bc.left);

    switch left_bc
        case 'non-slip'
            % normal velocity u = 0 at wall
            u(:, start) = 0.0;
            for k = 1:g
                u(:, start-k) = -u(:, start+k);
                v(:, start-k) = -v(:, start+k-1);
            end

        case 'free-slip'
            % normal velocity u = 0, tangential derivative dv/dn = 0
            u(:, start) = 0.0;
            for k = 1:g
                u(:, start-k) = -u(:, start+k);
                v(:, start-k) =  v(:, start+k-1);
            end

        otherwise
            error('Unsupported left velocity BC: %s', bc.left);
    end

    % ---------------- right boundary ----------------
    right_bc = normalize_bc(bc.right);

    switch right_bc
        case 'non-slip'
            u(:, endx+1) = 0.0;
            for k = 1:g
                if endx+1+k <= Nx_u
                    u(:, endx+1+k) = -u(:, endx+1-k);
                end
                if endx+k <= Nx_v
                    v(:, endx+k) = -v(:, endx-k+1);
                end
            end

        case 'free-slip'
            u(:, endx+1) = 0.0;
            for k = 1:g
                if endx+1+k <= Nx_u
                    u(:, endx+1+k) = -u(:, endx+1-k);
                end
                if endx+k <= Nx_v
                    v(:, endx+k) =  v(:, endx-k+1);
                end
            end

        otherwise
            error('Unsupported right velocity BC: %s', bc.right);
    end
end

%% ---------------- y-direction boundaries: down/up ----------------
if is_periodic_pair(bc.down, bc.up)

    % u-face periodic wrap
    for k = 1:g
        u(start-k, :) = u(endy-k+1, :);
        u(endy+k, :)  = u(start+k-1, :);
    end

    % v-face periodic wrap
    for k = 1:g
        v(start-k, :)   = v(endy+1-k, :);
        v(endy+1+k, :) = v(start+k, :);
    end

else
    % ---------------- down boundary ----------------
    down_bc = normalize_bc(bc.down);

    switch down_bc
        case 'non-slip'
            v(start, :) = 0.0;
            for k = 1:g
                u(start-k, :) = -u(start+k-1, :);
                v(start-k, :) = -v(start+k, :);
            end

        case 'free-slip'
            v(start, :) = 0.0;
            for k = 1:g
                u(start-k, :) =  u(start+k-1, :);
                v(start-k, :) = -v(start+k, :);
            end

        otherwise
            error('Unsupported down velocity BC: %s', bc.down);
    end

    % ---------------- up boundary ----------------
    up_bc = normalize_bc(bc.up);

    switch up_bc
        case 'non-slip'
            v(endy+1, :) = 0.0;
            for k = 1:g
                if endy+k <= Ny_u
                    u(endy+k, :) = -u(endy-k+1, :);
                end
                if endy+1+k <= Ny_v
                    v(endy+1+k, :) = -v(endy+1-k, :);
                end
            end

        case 'free-slip'
            v(endy+1, :) = 0.0;
            for k = 1:g
                if endy+k <= Ny_u
                    u(endy+k, :) =  u(endy-k+1, :);
                end
                if endy+1+k <= Ny_v
                    v(endy+1+k, :) = -v(endy+1-k, :);
                end
            end

        otherwise
            error('Unsupported up velocity BC: %s', bc.up);
    end
end

end

%% ========================================================================
% Helpers
%% ========================================================================
function bc = get_bc(grid)

bc = grid.bc;

bc.left  = normalize_bc(bc.left);
bc.right = normalize_bc(bc.right);
bc.down  = normalize_bc(bc.down);
bc.up    = normalize_bc(bc.up);

if xor(strcmp(bc.left,'periodic'), strcmp(bc.right,'periodic'))
    error('Periodic BC must be paired: left and right should both be periodic.');
end

if xor(strcmp(bc.down,'periodic'), strcmp(bc.up,'periodic'))
    error('Periodic BC must be paired: down and up should both be periodic.');
end

end

function s = normalize_bc(s)

s = lower(strtrim(s));

switch s
    case 'non-slip'
        s = 'non-slip';

    case 'free-slip'
        s = 'free-slip';

    case 'periodic'
        s = 'periodic';

    otherwise
        error('Unknown boundary condition type: %s', s);
end

end

function tf = is_periodic_pair(a, b)
tf = strcmp(a, 'periodic') && strcmp(b, 'periodic');
end

function A = fill_scalar_corners(A, start, endx, endy, Ny, Nx)

% 简单角点处理：复制最近有效角点
A(1:start-1, 1:start-1)       = A(start, start);
A(1:start-1, endx+1:Nx)       = A(start, endx);
A(endy+1:Ny, 1:start-1)       = A(endy, start);
A(endy+1:Ny, endx+1:Nx)       = A(endy, endx);

end