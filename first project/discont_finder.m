function [J0_discontinuities] = discont_finder(xIn1, xIn2, m)

    J0_discontinuities = zeros(1, m-10000);
    flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    epsilon = 0.001;

    tic;
    for i = 1:m-10000
        while(flag == 0)
            if counter > 2000
                error("Too many iterations, we have a problem")
            end
            
            f_left = besselj(0,lower);
            
            x_new = (lower + upper)/2;
            f_xnew = besselj(0, x_new);
            
            if abs(f_xnew) < epsilon
                J0_discontinuities(i) = x_new;
                flag = 1;
            elseif f_left*f_xnew < 0
                upper = x_new;
            else
                lower = x_new;
            end
            
            counter = counter + 1;
            
        end
        
        flag = 0;
        counter = 0;
        
        lower = x_new + 0.1;
        upper = x_new + 0.1 + pi;
    end
    toc;
end