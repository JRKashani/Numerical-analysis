function [p] = vectorized_laplacian_mat(p)
%"Takes the grid dimensions and averaging coefficient as input parameters.
% Returns a sparse, square matrix ($N \times N$) representing the system
% of linear equations required to solve for the temperature distribution,
% handling both interior logic and boundary constraints simultaneously.

    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    denom = p.n_avrg; %if wanted to weight different values to different
    % sides
    K = (1 : n_rows * n_cols)';
    %a single vector repesenting each coardinate of the temperature matrix.
    
    col_vect = ceil(K / n_rows);
    row_vect = mod(K-1, n_rows) + 1;
    %matlab using k = (j-1)*n_rows + 1, this will extract the wanted i and
    %j
    
    %boolean vectors meant to allow manipulation of the main vector
    top          = (row_vect == 1);
    bottom       = (row_vect == n_rows);
    left_most    = (col_vect == 1);
    right_most   = (col_vect == n_cols);
    p.the_boundary = (top | bottom | left_most | right_most);
    
    
    %the sparse command works with a vector of rows' coordinates, and
    %another of columbs', and a third of values. the next commands will
    %built those vectors.
    
    main_diag = ones(n_rows*n_cols, 1) * denom;
    main_diag(p.the_boundary) = 1;
    
    inner_square = K(~p.the_boundary);
    side_points  = ones(length(inner_square), 1) * (-1);
        
    sp_value_vect = [main_diag; repmat(side_points, 4, 1)];

    
    sp_row_vect = [K; repmat(inner_square, 4, 1)];
    sp_col_vect = [K; inner_square - 1; inner_square + 1;...
                   inner_square - n_rows; inner_square + n_rows];

    
    p.laplacian_mat = sparse(sp_row_vect, sp_col_vect, sp_value_vect,...
                      n_rows*n_cols, n_rows*n_cols);
end

% the_boundary = (row_vect == 1) | (row_vect == n_rows) |...
%         (col_vect == 1) | (col_vect == n_cols);
%is_inner_square = ~the_boundary;
% above_inner_row = inner_square;
    % above_inner_col = inner_square - 1;
    % 
    % below_inner_row = inner_square;
    % below_inner_col = inner_square + 1;
    % 
    % left_side_row = inner_square;
    % left_side_col = inner_square - n_rows;
    % 
    % right_side_row = inner_square;
    % right_side_col = inner_square + n_rows;
    %    
    % sp_row_vect = [main_diag_row; above_inner_row; below_inner_row;...
    %     left_side_row; right_side_row];
    % sp_col_vect = [main_diag_col; above_inner_col; below_inner_col;...
    %     left_side_col; right_side_col];
     % sp_value_vect = [main_diag; side_points; side_points; side_points;...
     %     side_points];
     % sp_row_vect = [K; inner_square; inner_square; inner_square;...
    %     inner_square];
     % main_diag_row = K;
    % main_diag_col = K;