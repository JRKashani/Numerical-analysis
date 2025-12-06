function [root_value] = project1(Bi, m)

    if ~isnumeric (Bi) || ~isscalar(Bi) || Bi >= 60 || Bi <= 0.1
        error("Your Biot number is irrelevant!");
    elseif ~isinteger(m) || ~isscalar(Bi) || m < 1
        error("The root number must be a natural number.");
    end

    epsilon = 10^(-6);
    flag = 0;
    counter = 0;

    if m = 1;
        root_value = Bi_section_root_finder(0.3, 2.4, Bi, epsilon);
    elseif m <= 1000
        J0_roots = readmatrix("roots_short.txt");
        higher_guess = J0_roots(m)-0.1;
        lower_guess = J0_roots(m-1)+0.1;
        root_value = ...
            Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    else
        J0_roots = readmatrix("roots_long.txt");
    end

    if m < 10000
        higher_guess = J0_roots(m)-0.1;
        lower_guess = J0_roots(m-1)+0.1;
        root_value = ...
            Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    else
        highest_known_discontinuity = J0_roots(10000);
        high_roots_array = zeros(1, m-10000);

        while(flag == 0)
            counter = counter + 1;
            if counter > 10000
                error("The main function is running too many iterrations");
            end

            temp_epsilon = 10^(-2);
            

        end
    end

    

    






end