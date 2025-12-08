function [root_value] = String_root_finder(Biot_number,...
    root_count, x_prev, x_n, epsilon, too_many_loops_text)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    flag = 0;
    counter = 0;
    

    while flag == 0
        counter = counter + 1;
        x_new = (f(x_n,Biot_number)*x_prev - f(x_prev, Biot_number)*x_n)...
            /(f(x_n, Biot_number)-f(x_prev, Biot_number));
        if counter > 10000
            error(too_many_loops_text);        
        elseif abs(f(x_new, Biot_number)) < epsilon
            root_value = x_new;
            flag = 1;
        else
            x_prev = x_n;
            x_n = x_new;
        end
    end
end