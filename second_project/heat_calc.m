function [temp_mat] = heat_calc(t_low,t_high)

    p.highest_temp = 100;
    p.lowest_temp = -30;
    p.mat_size_1 = 4;
    p.mat_size_2 = p.mat_size_1;
    p.n_parameters = p.mat_size_1 * p.mat_size_2;
    p.left_temp = 0;
    p.right_temp = -10;
    p.n_avrg = 4;
    
    validateattributes([t_low t_high], {'numeric'}, {'>=', ...
        p.lowest_temp,'<=', p.highest_temp, 'strictlyincreasing'});
    validateattributes([p.mat_size_1 p.mat_size_2], {'numeric'}, {'>=', ...
        2,'<=', 100});
    
    p.t_low = t_low;
    p.t_high = t_high;

    p.temp_mat = initial_mat_guess(p);

    p.laplacian_mat = initial_para_mat(p);





end