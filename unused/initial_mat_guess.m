function [Initialized_mat] = initial_mat_guess(p)
    
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;

    Initialized_mat = zeros(n_rows, n_cols);

    Initialized_mat(end, :) = p.t_low;
    Initialized_mat(1, :) = p.t_high;
    Initialized_mat(:, 1) = p.left_temp;
    Initialized_mat(:, end) = p.right_temp;

    y = linspace(0,1,n_rows).';
    x = linspace(0,1,n_cols);

    ii = 2:n_rows-1;
    jj = 2:n_cols-1;

     center = ...
     (1-y(ii))   .* Initialized_mat(1, jj)...
     + y(ii)     .* Initialized_mat(end, jj)...
     + (1-x(jj)) .* Initialized_mat(ii, 1)...
     + x(jj)     .* Initialized_mat(ii, end);
    
     Initialized_mat(ii, jj) = center / 2;
end

    % for i = 2:n_rows-1
    %     for j = 2:n_cols-1
    %         Initialized_mat(i,j) = ...
    %             ((i*Initialized_mat(end, j) + (n_rows - i)*...
    %             Initialized_mat(1, j))/(n_rows-1))/2 + ...
    %             ((j*Initialized_mat(i, end) + ...
    %             (n_cols - j)*Initialized_mat(i,1))/(n_cols-1))/2;
    %     end
    % end

        % vertical_cont = (i .* Initialized_mat(end, j) + (n_rows - i) .*...
    %     Initialized_mat(1, j))/(n_rows-1);
    % horizontal_cont = ((j .* Initialized_mat(i, end) +...
    %     (n_cols - j) .* Initialized_mat(i,1))/(n_cols-1));

       % Initialized_mat(i, j) = (vertical_cont + horizontal_cont) / 2;