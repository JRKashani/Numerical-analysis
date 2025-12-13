function benchmark_sparse_creation()
    % Define grid sizes to test
    sizes = [100, 500, 2000]; % 2000x2000 = 4 million nodes
    
    fprintf('------------------------------------------------------------\n');
    fprintf('| Grid Size (NxN) | Nodes (N) | Triplet Time (s) | spdiags Time (s) | Speedup |\n');
    fprintf('------------------------------------------------------------\n');

    for n = sizes
        p.mat_size_1 = n;
        p.mat_size_2 = n;
        p.n_avrg = 4;
        
        % Measure Triplet Method (Your Vectorized Code)
        func_triplet = @() create_triplets(p);
        time_triplet = timeit(func_triplet);
        
        % Measure spdiags Method
        func_spdiags = @() create_spdiags(p);
        time_spdiags = timeit(func_spdiags);
        
        % Calculate Speedup
        speedup = time_triplet / time_spdiags;
        
        fprintf('| %4d x %-4d   | %9d | %14.5f   | %14.5f   | %6.2fx |\n', ...
            n, n, n*n, time_triplet, time_spdiags, speedup);
    end
    fprintf('------------------------------------------------------------\n');
end

%% --- Method 1: Your Triplet Approach (Optimized) ---
function para_mat = create_triplets(p)
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    denom = p.n_avrg;
    N = n_rows * n_cols;
    
    K = (1:N)';
    col_vect = ceil(K / n_rows);
    row_vect = mod(K-1, n_rows) + 1;
    is_boundary = (row_vect == 1) | (row_vect == n_rows) | ...
                  (col_vect == 1) | (col_vect == n_cols);
    
    % Main diagonal
    main_diag_vals = ones(N, 1) * denom;
    main_diag_vals(is_boundary) = 1;
    
    % Neighbors
    inner_indices = K(~is_boundary);
    num_inner = length(inner_indices);
    side_points = ones(num_inner, 1) * -1;
    
    sp_row = [K; repmat(inner_indices, 4, 1)];
    sp_col = [K; inner_indices-1; inner_indices+1; ...
              inner_indices-n_rows; inner_indices+n_rows];
    sp_val = [main_diag_vals; repmat(side_points, 4, 1)];
    
    para_mat = sparse(sp_row, sp_col, sp_val, N, N);
end

%% --- Method 2: spdiags Approach ---
function para_mat = create_spdiags(p)
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    denom = p.n_avrg;
    N = n_rows * n_cols;
    
    K = (1:N)';
    col_vect = ceil(K / n_rows);
    row_vect = mod(K-1, n_rows) + 1;
    is_boundary = (row_vect == 1) | (row_vect == n_rows) | ...
                  (col_vect == 1) | (col_vect == n_cols);
    
    % Main Diagonal
    main_diag = ones(N, 1) * denom;
    main_diag(is_boundary) = 1;
    
    % Off-Diagonals (Initialize to -1)
    off_diag = ones(N, 1) * -1;
    
    % Mask off-diagonals: If a node is a boundary, it has 0 connection to others
    off_diag(is_boundary) = 0; 
    
    % Furthermore, spdiags wraps columns. We must break "false" connections
    % However, since boundary nodes are forced to Identity (1 on diag, 0 off),
    % the wrapping connections are naturally zeroed out by the mask above.
    
    % Diagonals: [Left (-nrows), Up (-1), Main (0), Down (+1), Right (+nrows)]
    B = [off_diag, off_diag, main_diag, off_diag, off_diag];
    d = [-n_rows, -1, 0, 1, n_rows];
    
    para_mat = spdiags(B, d, N, N);
end