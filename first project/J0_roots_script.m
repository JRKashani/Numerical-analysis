lower = 0.1;
upper = pi + 0.2;
counter = 0;
flag = 0;
epsilon = 1e-6;

J0_roots_short = zeros(1, 1000);
J0_roots_long = zeros(1, 10000);

for i = 1:1000
    while(flag == 0)
        if counter > 2000
            error("Too many iterations, we have a problem")
        end

        f_left = besselj(0,lower);

        x_new = (lower + upper)/2;
        f_xnew = besselj(0, x_new);

        if abs(f_xnew) < epsilon
            J0_roots_short(i) = x_new;
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

writematrix(J0_roots_short, 'roots_short.txt');

lower = 0.1;
upper = pi + 0.1;

for i = 1:10000
    while(flag == 0)
        if counter > 2000
            error("Too many iterations, we have a problem")
        end
        
        f_left = besselj(0,lower);
        
        x_new = (lower + upper)/2;
        f_xnew = besselj(0, x_new);
        
        if abs(f_xnew) < epsilon
            J0_roots_long(i) = x_new;
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

writematrix(J0_roots_long, 'roots_long.txt');

