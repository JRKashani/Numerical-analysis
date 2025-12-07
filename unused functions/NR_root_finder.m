function [root_value] = NR_root_finder(xIn, Bi_num)

    tic;
    flag = 0;
    counter = 0;
    x_n = xIn;
    epsilon = 10^(-5);
        
    while(flag == 0)
        if counter > 10000
            error("Too many iterations, we have a problem")
        end
        
        f_x = f(x_n, Bi_num);
        df_x = derv_f(x_n, Bi_num); 
        
        x_new = x_n - f_x/df_x;
        
        f_xnew = f(x_new, Bi_num);
        
        if abs(f_xnew) < epsilon
            root_value = x_new;
            flag = 1;
        else
            x_n = x_new;
        end
        disp([counter, x_new])
        counter = counter + 1;        
    end
    toc;
end

