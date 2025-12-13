
function benchmark_laplacian_builders(p)

    % Problem sizes to test (feel free to adjust)
    sizes = [50, 100, 200, 400];   % grid is size x size

    fprintf("Comparing performance of laplacian builders:\n");
    fprintf("----------------------------------------------------------\n\n");

    for s = sizes

        fprintf("Grid: %d x %d  (N = %d)\n", s, s, s*s);

        % Build input structure p
        p.mat_size_1 = s;
        p.mat_size_2 = s;
        p.n_avrg = 4;   % standard Laplacian weight

        % --- Vectorized builder ---
        t1 = tic;
        A = vectorized_laplacian_mat(p);
        t_vec = toc(t1);

        % --- spdiags builder (corrected) ---
        t2 = tic;
        B = laplacian_mat_diagonalised(p);
        t_spdiag = toc(t2);

        % --- Check correctness ---
        equal_structure = isequal(A ~= 0, B ~= 0);
        equal_values = norm(A - B, 'fro');

        fprintf("Vectorized build time:   %.6f sec\n", t_vec);
        fprintf("spdiags build time:      %.6f sec\n", t_spdiag);
        fprintf("Sparsity pattern match:  %s\n", logical_to_str(equal_structure));
        fprintf("Matrix difference norm:  %.3e\n", equal_values);

        fprintf("----------------------------------------------------------\n\n");
    end
end

function s = logical_to_str(x)
    if x
        s = "YES";
    else
        s = "NO";
    end
end
