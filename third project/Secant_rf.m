function [right_ddf_at_zero, counter] = Secant_rf(initial_guess)
    %The function utilize the secant root finding method to find the
    %initial value of f''(0) that will get f'(inf) = 1 as required

    func = @Blasius_Wrapping;
    
    %initializing the parameters
    flag = 0;
    upper_limit = 1;
    lower_limit = 0;
    epsilon = 10^-5;
    counter = 0;    
    x_old = 0;
    if initial_guess <= lower_limit || initial_guess > upper_limit
        x_curr = rand();
    else
        x_curr = initial_guess;
    end
    f1 = func(x_old);
    f2 = func(x_curr);
    f3 = 0;
    x_new = 0;

    while flag == 0
        counter = counter + 1;

        if counter > 1000
            error("Too many iterations!");
        end        
        
        %the secant formula
        x_new = x_curr -  f2*(x_curr - x_old)/(f2 - f1);

        f3 = func(x_new);

        if abs(f3) < epsilon
            flag = 1;
            right_ddf_at_zero = x_new;
        else
            x_old = x_curr;
            x_curr = x_new;
            f1 = f2;
            f2 = f3;
        end
    end
end