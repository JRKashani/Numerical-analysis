function [root_value,counter] = NR_root_finder(Biot_number,root_count, ...
    first_guess, epsilon)
%the function will find the root value of an accesory term for
% transient heat heat conduction of cylindrical items.

%   The function is guessing a number, and using the Newton-Rapson
%   method to find the nearest root of the function.

%initilizers
    too_many_loops_text = ['We are passed the 10000th loop, it seems'...
        'we have a problem'];
    counter = 0;
    flag = 0;
    x = first_guess;
    func_value = x*(besselj(1, x)/besselj(0, x)) - Biot_number;
        
    %%epsilon = 1e-3; % Define a small tolerance for convergence
    while flag == 0
        if abs(func_value) < epsilon
          root_value = x;
          break;
        elseif counter > 10000
           error(too_many_loops_text);
        end
        %in case of an unresoled question, the counter will prevent from
        %infinte loop
        counter = counter + 1;
        %getting the bessel functions values into compact form, for
        %improved readabillty
        J0 = besselj(0, x);
        J1 = besselj(1, x);
        J2 = besselj(2, x);
        %calculationg the function and derived values:
        func_value = x*(J1/J0) - Biot_number;
        deriv_value = (J1/J0)*(1+(x*(J1/J0))) - x/2 - J2/(2*J0);
        %the Newton-Raphson method core
        y = x - (func_value/deriv_value);
        x = y;
        
    end
    disp(counter);
end