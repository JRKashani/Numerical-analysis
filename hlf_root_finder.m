function [root_value,counter] = hlf_root_finder(Biot_number, lower_guess...
    , higher_guess, epsilon, too_many_loops_text)
%the function will find the root value of an accesory term for
% transient heat heat conduction of cylindrical items.

%   The function is guessing a number, and using the halfening
%   method to find the nearest root of the function.

%initilizers
    %too_many_loops_text = ['We are passed the 10000th loop, it seems we '...
     %      'have a problem'];
    counter = 0;
    flag = 0;
    func_value=10;
    %epsilon = 1e-3; % Define a small tolerance for convergence
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
        new = (lower_guess + higher_guess)/2;
        f_lower = func_calc(lower_guess, Biot_number);
        f_new = func_calc(new, Biot_number);
        if f_new*f_lower < 0
            higher_guess = new;
        elseif f_new*f_lower == 0
            root_value = new;
            flag = 1;
        else
            lower_guess = new;
        end
        
    end
    disp(counter);
end