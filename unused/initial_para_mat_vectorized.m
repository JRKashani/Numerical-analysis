function [para_mat] = initial_para_mat_vectorized(p)
    rows = p.mat_size_1;
    cols = p.mat_size_2;
    N = rows * cols;
    denom = p.n_avrg;

    % --- Phase 1 & 2: Coordinates and Masks ---
    K = (1:N)'; % Master index col vector
    I_vec = mod(K - 1, rows) + 1;
    J_vec = ceil(K / rows);

    is_boundary = (I_vec == 1) | (I_vec == rows) | ...
                  (J_vec == 1) | (J_vec == cols);
    is_interior = ~is_boundary;

    % --- Phase 3: Building Triplet Chunks ---
    % 3A. Diagonal Connections (All points)
    diag_I = K;
    diag_J = K;
    diag_V = ones(N, 1) * denom;
    diag_V(is_boundary) = 1; % Overwrite boundaries using mask

    % 3B. Neighbor Connections (Interior only)
    K_in = K(is_interior); % Filtered index list
    num_int = length(K_in);
    neg_ones = ones(num_int, 1) * -1;

    % Define offsets for column-major indexing
    up_J    = K_in - 1;
    down_J  = K_in + 1;
    left_J  = K_in - rows;
    right_J = K_in + rows;

    % --- Phase 4: Assembly ---
    % Stack everything vertically.
    % Note: K_in is used repeatedly as the 'I' coordinate for neighbors
    final_I = [diag_I; K_in;   K_in;   K_in;   K_in];
    final_J = [diag_J; up_J;   down_J; left_J; right_J];
    final_V = [diag_V; neg_ones; neg_ones; neg_ones; neg_ones];

    para_mat = sparse(final_I, final_J, final_V, N, N);
end