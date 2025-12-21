function [inversed_mat] = inverse_this_diag_matrix(diag_mat_for_inversion)
%INVERSE_THIS_DIAG_MATRIX receives a diagonal matrix and returns its
%inverse, calculated element-wise without using built-in matrix inversion
% functions.
    
    % Initializing parameters
    n_rows = size(diag_mat_for_inversion, 1);
    aux_vect = zeros(n_rows, 1);

    % Iterating through diagonal elements
    for i = 1 : n_rows
        % Inverting each element scalar-wise (1/x) to build the
        % diagonal vector
        aux_vect(i) = 1 / diag_mat_for_inversion(i, i);
    end
    % Reconstructing the sparse diagonal matrix from the inverted vector
    inversed_mat = spdiags(aux_vect, 0, n_rows, n_rows);
end