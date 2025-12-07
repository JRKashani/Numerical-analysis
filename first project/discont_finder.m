function [x_prev, x_m] = discont_finder(xIn1, xIn2, m)
%This function will return the (m-1)-th and m-th roots of the
%zeroeth bessel function after recieving 2 initial values to find the
%first root between them, and the wanted m.
    
    %initilizing the parameters, and aux variables.
    flag = 0;
    sec_flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    x_prev = 0;
    x_m = 0;
    epsilon = 0.0001;
    if m <= 10000
        m = m + 10000;
    end

    for i = 1:m-10000
        %the loop will run m-10000 times to return the m-th root while
        %constantly keeping the previous root.
        while(flag == 0)
            if counter > 2000
                error("Too many iterations, we have a problem")
                %preventing infinte loop, the counter is adding in the end
                %of it
            end

            if i > (m - 10000 - 2)
                epsilon = 10^(-7);
                %improving the accuracy in the last roots, the ones that
                %will be used.
            end
            
            if sec_flag == 0
                f_lower = besselj(0,lower);
                %saving the new evalution of the function each time the
                %lower boundry isn't changing.
            end
            sec_flag = 0;    
            x_new = (lower + upper)/2;
            f_xnew = besselj(0, x_new);
            %calculating the value of the new x, assesing its function
            %value.
            
            if abs(f_xnew) < epsilon
                if i > (m - 10000 - 2)
                    x_prev = x_m;
                    x_m = x_new;
                end
                flag = 1;
                %each root being identified, but only the relevant one
                %stored
                
            elseif f_lower*f_xnew < 0
                upper = x_new;
                sec_flag = 1;
                %if not close enough to a root, cutting by half the range
                %and checking for different signs of the function at the
                %boundries.
            else
                lower = x_new;
            end
            
            counter = counter + 1;
            
        end

        %initilizing again for the next root
        
        flag = 0;
        sec_flag = 0;
        counter = 0;
        
        lower = x_new + 0.1;
        upper = x_new + 0.1 + pi;
    end
end