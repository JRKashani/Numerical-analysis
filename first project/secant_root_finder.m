function [root_value] = secant_root_finder(Inx1, Inx2, Bi_num, epsilon)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    flag = 0;
    counter = 0;
    f1 = 0;
    f2 = 0;
    f3 = 0;

    x_prev = Inx1;
    x_n = Inx2;


    while(flag == 0)
        if counter > 10000
            error("Too many iterations, we have a problem")
        end

        f1 = f(x_prev, Bi_num);
        f2 = f(x_n, Bi_num);

        x_new = x_n -  f2*(x_n - x_prev)/(f2 - f1);

        f3 = f(x_new, Bi_num);

        if abs(f3) < epsilon
            flag = 1;
            root_value = x_new;
        else
            x_prev = x_n;
            x_n = x_new;
        end

        counter = counter + 1;
    end

end

