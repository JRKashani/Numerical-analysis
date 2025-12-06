function [root_value] = project1(Bi, m)

    if ~isnumeric (Bi) || ~isscalar(Bi) || Bi >= 60 || Bi <= 0.1
        error("Your Biot number is irrelevant!");
    elseif ~isnumeric(m) || m <= 0 || m ~= floor(m)
        error("The root number must be a natural number.");
    end

    epsilon = 10^(-6);

    if m == 1
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

    if m > 1000 && m <= 10000
        higher_guess = J0_roots(m)-0.1;
        lower_guess = J0_roots(m-1)+0.1;
        root_value = ...
            Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    elseif m > 10000
        J0_discontinuity = discont_finder(J0_roots(10000) + 0.01, ...
            J0_roots(10000) + 0.01 + pi, m);

        lower_discon = J0_discontinuity(m - 10000 - 1) + 0.1;
        higher_discon = J0_discontinuity(m - 10000) - 0.1;
        root_value = Bi_section_root_finder(lower_discon, higher_discon,...
            Bi, epsilon);
    end
end