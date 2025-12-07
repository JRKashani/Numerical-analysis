function [root_value] = project1(Bi, m)
%The function will recive a Biot number and a natural number "m" and return
%the m-th root of the transient heat conduct equation of a cylindrical
%item.

%checking the input arguments for non-numeric or not in range user-input.
    if ~isnumeric (Bi) || ~isscalar(Bi) || Bi > 60 || Bi < 0.1
        error("Your Biot number is irrelevant!");
    elseif ~isnumeric(m) || m <= 0 || m ~= floor(m)
        error("The root number must be a natural number.");
    end

%defining the accuracy which we aim for.
    epsilon = 10^(-6);

%to save run time, we made 2 txt files containing the discontinuites of 
%the function, which are the roots of the zeroeth bessel function, more 
%specifically the first 100 roots, and the first 10000, and to find the
%m-th root of our function, we accses the m-th and (m-1)-th 
%discontinuity and search for a root in those boundaries, assuring we 
%won't be looking around a discontinuity in the function, and won't
%need to check "m" roots of the function.
    if m == 1
        root_value = Bi_section_root_finder(0.3, 2.4, Bi, epsilon);
    elseif m <= 100
        J0_roots = readmatrix("roots_short.txt");
        %given the Biot numbers, a distance of 0.1 from the discontinuty
        % will suffice for every root. 
        higher_guess = J0_roots(m) - 0.1; 
        lower_guess = J0_roots(m-1) + 0.1;
        root_value = ...
            Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    else
        J0_roots = readmatrix("roots_long.txt");
        %the readmatrix comand is slow for this size of a file, so we
        %preffered to split the roots files and open only the neccesary
        %one.
    end

    if m > 100 && m <= 10000
        higher_guess = J0_roots(m) - 0.1;
        lower_guess = J0_roots(m-1) + 0.1;
        root_value = ...
            Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    elseif m > 10000
        %if there is no other choise, we will continue mapping the
        %discontinuiues and find the root between the (m-1)-th and m-th
        %discontinuity points. this is obviuosly the longest and least
        %desired path to take. the "pi" distance is aprroximatly the
        %bessel function frequncy.

        [lower_discon, higher_discon] = discont_finder(J0_roots(10000)...
            + 0.01, J0_roots(10000) + 0.01 + pi, m);
        root_value = Bi_section_root_finder(lower_discon + 0.1,...
            higher_discon - 0.1, Bi, epsilon);
    end
    
end