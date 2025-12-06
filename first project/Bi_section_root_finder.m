function [root_value] = Bi_section_root_finder(xIn1, xIn2, Bi_num, epsilon)

    flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    root_value = 0;
            
    while(flag == 0)
        if counter > 10000
            error("Too many iterations, we have a problem")
        end
        
        f_left = f(lower, Bi_num);
        
        x_new = (lower + upper)/2;
        f_xnew = f(x_new, Bi_num);
        
        if abs(f_xnew) < epsilon
            root_value = x_new;
            flag = 1;
        elseif f_left*f_xnew < 0
            upper = x_new;
        else
            lower = x_new;
        end
        
        counter = counter + 1;
    end
end

