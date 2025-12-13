function [laplacian_mat] = slow_laplacian_mat(p)
    
    laplacian_mat = zeros(p.n_parameters, p.n_parameters);

    for i = 1:p.n_parameters
        if i <= p.mat_size_1 ||...
                i >= p.n_parameters - p.mat_size_1 ||...
                mod(i, p.mat_size_2) == 0 ||...
                mod(i, p.mat_size_2) == 1
            laplacian_mat(i,i) = 1;
        else
            laplacian_mat(i,i) = p.n_avrg;
            laplacian_mat(i,i - 1) = -1;
            laplacian_mat(i,i + 1) = -1;
            laplacian_mat(i,i + p.mat_size_1) = -1;
            laplacian_mat(i,i - p.mat_size_1) = -1;
        end        
    end

end