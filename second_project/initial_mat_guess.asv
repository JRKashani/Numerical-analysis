function [Initialized_mat] = initial_mat_guess(p)
    
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;

    Initialized_mat = zeros(n_rows, n_cols);

    Initialized_mat(end, :) = p.t_low;
    Initialized_mat(1, :) = p.t_high;
    Initialized_mat(:, 1) = p.left_temp;
    Initialized_mat(:, end) = p.right_temp;

    i = (2:n_rows - 1);
    j = (2:n_cols - 1);

    for i = 2:n_rows-1
        for j = 2:n_cols-1
            Initialized_mat(i,j) = ...
                ((i*Initialized_mat(end, j) + (n_rows - i)*...
                Initialized_mat(1, j))/(n_rows-1))/2 + ...
                ((j*Initialized_mat(i, end) + ...
                (n_cols - j)*Initialized_mat(i,1))/(n_cols-1))/2;
        end
    end

end