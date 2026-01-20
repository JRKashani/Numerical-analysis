function [p] = laplacian_solver_GS_4step_vectorized(p)
% LAPLACIAN_SOLVER_GS_4STEP_VECTORIZED solves the system using a 4-stage 
% update with pre-sliced matrices for maximum vectorization speed.

    % Unpacking
    linear_equations_matrix = p.laplacian_mat;
    x_curr = p.temp_mat(:);
    
    % Construct the RHS vector 'b'
    b_full = x_curr .* p.the_boundary(:);
    
    % Decomposition
    % D = Diagonal matrix
    % LU = Off-diagonal matrix (L + U)
    n = p.n_parameters;
    D_diag = diag(linear_equations_matrix);
    D = spdiags(D_diag, 0, n, n);
    LU_full = linear_equations_matrix - D; 
    D1_full = p.inv_D; 

    % --- PRE-CALCULATION / PRE-SLICING ---
    % We pre-allocate cells to hold the 4 distinct parts of the system.
    % This moves the expensive indexing operation out of the iterative loop.
    indices  = cell(4, 1);
    LU_parts = cell(4, 1);
    b_parts  = cell(4, 1);
    D1_parts = cell(4, 1);
    
    all_indices = (1:n)';
    
    for k = 1:4
        % Define the group indices (0, 1, 2, 3 mod 4)
        % We use k-1 to match the mod 4 result (0..3)
        idx = all_indices(mod(all_indices, 4) == (k-1));
        
        indices{k}  = idx;
        
        % CRITICAL STEP: Slice the sparse matrix ONCE here.
        % Inside the loop, we will only do multiplication.
        LU_parts{k} = LU_full(idx, :); 
        
        b_parts{k}  = b_full(idx);
        
        % Slice the inverse diagonal for this group
        D1_parts{k} = D1_full(idx, idx);
    end
    
    % Initialize loop variables
    counter = 0;
    flag = 0;
    x_prev_check = x_curr; % Snapshot for convergence check
    
    while flag == 0
        counter = counter + 1;
        
        % --- 4-STAGE VECTORIZED UPDATE ---
        for k = 1:4
            % 1. Calculate row sums for the current group using the 
            %    CURRENT full x vector (which includes updates from k-1).
            %    This is a fast Sparse-Matrix * Vector operation.
            row_sum = LU_parts{k} * x_curr;
            
            % 2. Update only the relevant indices in x_curr
            x_curr(indices{k}) = D1_parts{k} * (b_parts{k} - row_sum);
        end

        % --- CONVERGENCE CHECK ---
        % Checked every 10 iterations to amortize cost
        if mod(counter, 10) == 0
            if counter > 100000
                error("exceeded the maximum number of iterations");
            end
            
            % Calculate diff against the snapshot from 10 iterations ago
            diff_vect = abs(x_curr - x_prev_check);
            
            if max(diff_vect) < p.epsilon
                flag = 1;
            elseif max(diff_vect) > 100
                error("starting to diverge, no point of continuing");
            end
            
            % Update snapshot
            x_prev_check = x_curr;
        end
    end
    
    p.temp_mat = reshape(x_curr, p.mat_size_1, p.mat_size_2);
end