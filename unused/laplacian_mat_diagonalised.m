% function [laplacian_mat] = laplacian_mat_diagonalised(p)
% 
%     n_rows = p.mat_size_1;
%     n_cols = p.mat_size_2;
%     N = n_rows*n_cols;
%     denom = p.n_avrg; %if wanted to weight different values to different
%     % sides
%     K = (1:N)';
%     %a single vector repesenting each coardinate of the temperature matrix.
% 
%     col_vect = ceil(K / n_rows);
%     row_vect = mod(K-1, n_rows) + 1;
%     %matlab using k = (j-1)*n_rows + 1, this will extract the wanted i and
%     %j
% 
%     the_boundary = (row_vect == 1) | (row_vect == n_rows) |...
%                    (col_vect == 1) | (col_vect == n_cols);
% 
%     main_diag = ones(N, 1) * denom;
%     main_diag(the_boundary) = 1;
%     side_points = ones(N, 1) * (-1);
% 
%     side_points(the_boundary) = 0;
%     all_diags = [side_points, side_points, main_diag, side_points,...
%                  side_points];
%     diag_offset = [(-n_rows), (-1), 0, 1, n_rows];
% 
%     laplacian_mat = spdiags(all_diags, diag_offset, N, N);
% end
function L = laplacian_mat_diagonalised(p)
    % Builds exactly the same matrix as vectorized_laplacian_mat(p)
    % using spdiags, including repeated neighbor contributions.

    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    N = n_rows * n_cols;
    denom = p.n_avrg;

    K = (1:N)';

    col_vect = ceil(K / n_rows);
    row_vect = mod(K-1, n_rows) + 1;

    % boundary nodes
    the_boundary = (row_vect == 1) | (row_vect == n_rows) | ...
                   (col_vect == 1) | (col_vect == n_cols);

    % interior nodes
    interior = ~the_boundary;

    % main diagonal: 1 for boundary, denom + sum of off-diagonals for interior
    main_diag = ones(N,1) * denom;
    main_diag(the_boundary) = 1;

    % Off-diagonals: -1 for interior nodes
    sp_up   = zeros(N,1); sp_up(interior) = -1;
    sp_down = zeros(N,1); sp_down(interior) = -1;
    sp_left = zeros(N,1); sp_left(interior) = -1;
    sp_right= zeros(N,1); sp_right(interior) = -1;

    % Adjust main diagonal for interior nodes to account for 4 neighbors
    main_diag(interior) = main_diag(interior) + 4*(-1);

    % Build spdiags
    all_diags = [sp_up, sp_left, main_diag, sp_right, sp_down];
    diag_offset = [-n_rows, -1, 0, 1, n_rows];

    L = spdiags(all_diags, diag_offset, N, N);
end
