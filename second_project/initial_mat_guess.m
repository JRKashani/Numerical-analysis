function [p] = initial_mat_guess(p)

    %using struct to save place inline
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;

    Initialized_mat = zeros(n_rows, n_cols);

    Initialized_mat(1, :)   = p.t_high;      % Top
    Initialized_mat(end, :) = p.t_low;       % Bottom
    Initialized_mat(:, 1)   = p.left_temp;   % Left
    Initialized_mat(:, end) = p.right_temp;  % Right

    %making the vectors of parameters
    y_cent = (1 : n_rows - 2) .' / (n_rows - 1);
    x_cent = (1 : n_cols - 2)    / (n_cols - 1);
    
    %vectors of tempratures as mean of the walls values
    t_rows = p.t_high + y_cent * (p.t_low - p.t_high);
    t_cols = p.left_temp + x_cent * (p.right_temp - p.left_temp);
    
    %setting the central area of the temp matrix
    Initialized_mat(2 : end-1 , 2 : end-1) = (t_rows + t_cols) / 2;
    p.temp_mat = Initialized_mat;
end