function [inversed_mat] = inverse_this_diag_matrix(diag_mat_for_inversion)
    
    n_rows = size(diag_mat_for_inversion, 1);
    aux_vect = zeros(n_rows, 1);

    for i = 1 : n_rows
        aux_vect(i) = 1 / diag_mat_for_inversion(i, i);
    end

    inversed_mat = spdiags(aux_vect, 0, n_rows, n_rows);
end