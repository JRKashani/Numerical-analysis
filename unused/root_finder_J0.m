function [root] = root_finder_J0(lower, upper)

    counter = 0;
    flag = 0;
    epsilon = 1e-10;
    if besselj(0,lower)*besselj(0,upper) > 0
        error("wrong limits.")
    end

    
    while(flag == 0)
        counter = counter + 1;
        new_x = (upper+lower)/2;

        if counter > 10000
            error("Too many iterations!");
        elseif abs(besselj(0, new_x)) < epsilon
            root = new_x;
            flag = 1;
        elseif besselj(0,lower)*besselj(0,new_x) < 0
            upper = new_x;
        else
            lower = new_x;
        end
    end
    %disp(counter);

end