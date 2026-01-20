function [p] = vectorized_laplacian_mat(p)
% VECTORIZED_LAPLACIAN_MAT Constructs the discrete 2D Laplacian operator.
% The function builds a sparse (N x N) matrix corresponding to the finite-
% difference Laplacian on a rectangular grid, including enforcement of
% Dirichlet boundary conditions via identity rows.

    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    denom = p.n_avrg; %if wanted to weight different values to different
    % sides
    K = (1 : n_rows * n_cols)';
    % Linear index vector corresponding to each grid point in column-major
    % order
    
    col_vect = ceil(K / n_rows);
    row_vect = mod(K-1, n_rows) + 1;
    % MATLAB uses column-major ordering: k = (j-1)*n_rows + i.
    % The following expressions recover (i,j) grid indices from k.
    
    %boolean vectors meant to allow manipulation of the main vector
    top            = (row_vect == 1      );
    bottom         = (row_vect == n_rows );
    left_most      = (col_vect == 1      );
    right_most     = (col_vect == n_cols );
    p.the_boundary = (top | bottom | left_most | right_most);
    
    
    % The sparse constructor requires vectors of row indices, column 
    % indices, and corresponding nonzero values; the following code 
    % assembles these.

    
    % Main diagonal entries: 4 for interior nodes (Laplacian stencil),
    % and 1 for boundary nodes to enforce Dirichlet conditions via
    % identity rows
    main_diag = ones(n_rows*n_cols, 1) * denom;
    inv_diag  = ones(n_rows*n_cols, 1) * 1/denom;
    main_diag(p.the_boundary) = 1;
    inv_diag(p.the_boundary) = 1;
    
    % Off-diagonal entries corresponding to nearest-neighbor interactions
    % (−1 for each interior neighbor in the 5-point stencil)
    inner_square = K(~p.the_boundary);
    side_points  = ones(length(inner_square), 1) * (-1);
    
    % Assemble value vector for the sparse multi-diagonal matrix:
    % one main diagonal and four neighbor diagonals for interior nodes
    sp_value_vect = [main_diag; repmat(side_points, 4, 1)];

    
    % Row and column index vectors for sparse assembly.
    % Each interior node contributes four additional entries corresponding
    % to its left, right, top, and bottom neighbors.
    sp_row_vect = [K; repmat(inner_square, 4, 1)];
    sp_col_vect = [K; inner_square - 1; inner_square + 1;...
                   inner_square - n_rows; inner_square + n_rows];

    % Sparse storage significantly reduces memory usage and computation cost
    p.laplacian_mat = sparse(sp_row_vect, sp_col_vect, sp_value_vect,...
                      n_rows*n_cols, n_rows*n_cols);
    p.inv_D = sparse(K, K, inv_diag, n_rows*n_cols, n_rows*n_cols);
end

% --- Legacy / alternative sparse assembly approach (not used) ---

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
