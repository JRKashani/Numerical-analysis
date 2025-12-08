function [root_value] = Bi_section_root_finder(xIn1, xIn2, Bi_num, epsilon)
%This function utelize the Bi section method for finding a root of a
%function in a limited range - the function will recive boundries, biot
%number and tolerance, and return an 'x' which its function value is within
%an epsilon from zero.

    %initilizing the parameters and aux variables
    flag = 0;
    sec_flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    root_value = 0;

    if f(lower, Bi_num)*f(upper, Bi_num) > 0
        upper = lower + pi/2 + 0.1;
        %for cases the boundries weren't defined correctly, this will
        %provide a result, although less reliable one.
    end
            
    while(flag == 0)
        if counter > 10000
            error("Too many iterations, we have a problem")
            %prevents an endless loop
        end
        
        %to prevent it from running the besselj function each time, it will
        %do only if the lower x value had been changed.
        if sec_flag == 0
            f_lower = f(lower, Bi_num);
        end
        sec_flag = 0;
        
        x_new = (lower + upper)/2;
        %x_new will be the average of the prevoius 'x's.
        f_xnew = f(x_new, Bi_num);
        
        if abs(f_xnew) < epsilon
            %if the new x value isn't close enough, we will choose the x
            %that yield result with different sign than the current.
            root_value = x_new;
            flag = 1;
        elseif f_lower*f_xnew < 0
            upper = x_new;
            sec_flag = 1;
        else
            lower = x_new;
        end
        
        counter = counter + 1;
    end
end

