function [Initialized_mat] = initial_mat_guess(p)
    % Unpack dimensions
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    
    % 1. Create the base matrix with Boundary Conditions
    Initialized_mat = zeros(n_rows, n_cols);
    
    % Set Boundaries (Dirichlet)
    % Note: The order matters for corners. The last one written "wins".
    Initialized_mat(1, :)   = p.t_high;      % Top
    Initialized_mat(end, :) = p.t_low;       % Bottom
    Initialized_mat(:, 1)   = p.left_temp;   % Left
    Initialized_mat(:, end) = p.right_temp;  % Right
    
    % 2. Calculate Interior Guess (Vectorized)
    % We only calculate for the interior points (2 to N-1)
    
    % Normalized coordinates for the INTERIOR only
    % (Avoiding 0 and 1 because those are boundaries)
    y_inner = linspace(0, 1, n_rows).'; 
    y_inner = y_inner(2:end-1);         % Vertical vector (N-2 x 1)
    
    x_inner = linspace(0, 1, n_cols);
    x_inner = x_inner(2:end-1);         % Horizontal vector (1 x M-2)
    
    % Vertical Linear Interpolation (Top to Bottom)
    % Formula: (1-y)*Top + y*Bottom
    t_vert = (1 - y_inner) * p.t_high + y_inner * p.t_low;
    
    % Horizontal Linear Interpolation (Left to Right)
    % Formula: (1-x)*Left + x*Right
    t_horz = (1 - x_inner) * p.left_temp + x_inner * p.right_temp;
    
    % 3. Combine and Fill
    % Implicit expansion automatically creates the grid:
    % (N-1 x 1) + (1 x M-1) = (N-1 x M-1) Matrix
    Initialized_mat(2:end-1, 2:end-1) = (t_vert + t_horz) / 2;

end