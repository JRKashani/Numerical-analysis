function [temp_mat] = heat_calc(t_low,t_high)
% HEAT_CALC Solves the steady-state 2D heat equation on a square grid
% using Gauss–Seidel or Jacobi methods and returns the heat matrix.
%
%   temp_mat = HEAT_CALC(t_low, t_high)
%
% Inputs:
%   t_low   - Lower boundary temperature
%   t_high  - Upper boundary temperature
%
% Output:
%   temp_mat - heat matrix
    
    % Initialize solver parameter structure
    p.highest_temp = 100;
    p.lowest_temp = -30;
    p.mat_size_1 = 16;
    p.mat_size_2 = p.mat_size_1;
    p.n_parameters = p.mat_size_1 * p.mat_size_2;
    p.left_temp = 0;
    p.right_temp = -10;
    p.n_avrg = 4; % Number of neighboring points in 2D Laplacian stencil
    p.epsilon = 10 ^ (-6);
    
   % Validate input arguments and internal parameters
    validateattributes([t_low t_high], {'numeric'}, {'>=', ...
        p.lowest_temp,'<=', p.highest_temp, 'increasing'});
    validateattributes([p.mat_size_1 p.mat_size_2], {'numeric'}, {'>=', ...
        2,'<=', 100});
    if t_high == t_low
        error("Boundary temperatures must not be equal.");
    end
    
    %entering the validated data
    p.t_low = t_low;
    p.t_high = t_high;

    % To reduce convergence time, initialize the temperature field using
    % a weighted average of the boundary temperatures rather than zeros
    p = initial_mat_guess(p);
    
    % Construct the discrete Laplacian matrix (N x N), one row per grid
    % point. Boundary nodes are retained for indexing consistency; their
    % rows are replaced with identity rows to enforce Dirichlet boundary
    % conditions.
    p = vectorized_laplacian_mat(p);

    % Gauss–Seidel method solver.
    % Although GS converges in fewer iterations than Jacobi, for this
    % problem the matrix-based implementation is slower in MATLAB due to
    % reduced vectorization efficiency.
    %p = laplacian_solver_GS(p);

    %GS = p.temp_mat;
    
    % Jacobi method solver.
    % Due to its fully vectorized matrix formulation, Jacobi performs
    % better than GS in the MATLAB execution environment for this problem.

    p = laplacian_solver_Jac(p);

    Jac = p.temp_mat;

    temp_mat = Jac;
end