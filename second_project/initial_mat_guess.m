function [p] = initial_mat_guess(p)
% INITIAL_MAT_GUESS Generates an initial temperature field for the solver.
% The function accepts a parameter structure describing the domain and
% boundary conditions, and returns an initial guess that satisfies the
% boundaries and interpolates smoothly across the interior.

    % Extract grid dimensions from parameter structure
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;

    % Initialize temperature matrix and assign Dirichlet boundary values
    Initialized_mat = zeros(n_rows, n_cols);

    Initialized_mat(1, :)   = p.t_high;      % Top
    Initialized_mat(end, :) = p.t_low;       % Bottom
    Initialized_mat(:, 1)   = p.left_temp;   % Left
    Initialized_mat(:, end) = p.right_temp;  % Right

    % Construct normalized interior grid coordinates
    y_cent = (1 : n_rows - 2) .' / (n_rows - 1);
    x_cent = (1 : n_cols - 2)    / (n_cols - 1);
    
    % Linearly interpolated temperature profiles between opposing
    % boundaries
    t_rows = p.t_high + y_cent * (p.t_low - p.t_high);
    t_cols = p.left_temp + x_cent * (p.right_temp - p.left_temp);
    
    % This initialization provides a smooth starting field, reducing the
    % number of iterations required for convergence of iterative solvers

    % Initialize interior nodes as the average of horizontal and vertical
    % linear boundary interpolations
    Initialized_mat(2 : end-1 , 2 : end-1) = (t_rows + t_cols) / 2;
    p.temp_mat = Initialized_mat;
end