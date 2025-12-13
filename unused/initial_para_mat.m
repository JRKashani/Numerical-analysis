function [para_mat] = initial_para_mat(p)
%"Takes the grid dimensions and averaging coefficient as input parameters.
% Returns a sparse, square matrix ($N \times N$) representing the system
% of linear equations required to solve for the temperature distribution,
% handling both interior logic and boundary constraints simultaneously.
    
    %The struct is defined in the main function
    n_rows = p.mat_size_1;
    n_cols = p.mat_size_2;
    denom = p.n_avrg; %if wanted to weight different values to different
    % sides
    counter = 1;

    n_temp_points = n_rows * n_cols;
    nz_elements = (denom + 1) * (n_cols-1) * (n_rows-1) +...
        2 * (n_rows) + 2 * (n_cols - 2);

    vect_of_row_cord = ones(1, nz_elements);
    vect_of_col_cord = ones(1, nz_elements);
    vect_of_inserted_values = ones(1, nz_elements);

    for k = 1:n_temp_points
        % the index is calculated as k = n_rows*(j-1) + i so to extract j
        % and i, we will use modulo operator and ceil
        j = ceil(k / n_rows);
        i = mod(k-1, n_rows) + 1;

        top = (i == 1);
        bottom = (i == n_rows);
        left_most = (j == 1);
        right_most = (j == n_cols);

        if top || bottom || left_most || right_most
            vect_of_row_cord(counter) = k;
            vect_of_col_cord(counter) = k;
            vect_of_inserted_values(counter) = 1;
            counter = counter + 1;
            continue;
        else
            vect_of_row_cord(counter) = k;
            vect_of_col_cord(counter) = k;
            vect_of_inserted_values(counter) = denom;
            counter = counter + 1;
            
            sides = [k-1, k+1, k - n_rows, k + n_rows];

            for side = sides
                vect_of_row_cord(counter) = k;
                vect_of_col_cord(counter) = side;
                vect_of_inserted_values(counter) = -1;
                counter = counter + 1;
            end
        end
    end
    para_mat = sparse(vect_of_row_cord(1:counter-1),...
        vect_of_col_cord(1:counter-1),...
        vect_of_inserted_values(1:counter-1),...
        n_temp_points, n_temp_points);
end