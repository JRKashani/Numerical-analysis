function root_value = cordway(Bi,xIn1,xIn2)

    epsilon = 1e-4;
    x_prev = xIn1; %x_n-1
    x_current = xIn2; %x_n
    flag = 0;
    counter = 0;
    while flag == 0
        if counter > 1000000
            error('Too many iterations');
        else
            x_next = (f(x_current,Bi)*x_prev-f(x_prev,Bi)*x_current)/...
                (f(x_current,Bi)-f(x_prev,Bi));
        end
        if abs(f(x_next, Bi)) < epsilon
            flag = 1;
        else
            counter = counter + 1;
            x_prev = x_current;
            x_current = x_next;
        end
    end
    root_value = x_next;
end