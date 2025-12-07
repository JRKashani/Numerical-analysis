% RUN_PROJECT_COMPARISON
% Generates verification data and performance metrics for Project 1.
% Compares 'SolveBesselAI' (Hybrid Newton) vs. 'SolveBesselManual' (Bisection).

clc; clear; close all;

fprintf('==============================================================\n');
fprintf('               PROJECT 1: SOLVER COMPARISON                   \n');
fprintf('==============================================================\n');

% --- Configuration ---
% Define test cases: [Bi, m]
test_cases = [
    24.8, 1;   % Required Case 1
    24.8, 2;   % Required Case 2
    24.8, 5;   % Required Case 3
    0.5,  1;   % Small Bi
    60.0, 3;   % Large Bi
    10.0, 50;  % Large m (Stress Test)
];

fprintf('\n--- PART 1: ACCURACY & ROOT LOCATION ---\n');
fprintf('%-6s | %-4s | %-12s | %-12s | %-10s\n', ...
        'Bi', 'm', 'x (Manual)', 'x (AI)', 'Diff');
fprintf('--------------------------------------------------------------\n');

for i = 1:size(test_cases, 1)
    Bi = test_cases(i, 1);
    m  = test_cases(i, 2);
    
    % 1. Run AI Solver
    try
        t_start = tic;
        x_ai = SolveBesselAI(Bi, m);
        t_ai = toc(t_start);
    catch
        x_ai = NaN;
    end
    
    % 2. Run Manual Solver 
    try
        x_manual = project1(Bi, m);
    catch
        x_manual = NaN;
    end
    
    % 3. Calculate Difference
    diff_val = abs(x_ai - x_manual);
    
    fprintf('%6.1f | %4d | %12.6f | %12.6f | %10.2e\n', ...
            Bi, m, x_manual, x_ai, diff_val);
end

fprintf('\n\n--- PART 2: PERFORMANCE SPEED TEST (100,000 Runs) ---\n');
Bi_speed = 24.8; 
m_speed = 1;
num_runs = 100000;

fprintf('Benchmarking for Bi=%.1f, m=%d...\n', Bi_speed, m_speed);

% Measure Manual
tic;
for k = 1:num_runs
    dummy = project1(Bi_speed, m_speed);
end
time_manual = toc;

% Measure AI
tic;
for k = 1:num_runs
    dummy = SolveBesselAI(Bi_speed, m_speed);
end
time_ai = toc;

fprintf('Manual Code Total Time:  %8.4f seconds\n', time_manual);
fprintf('AI Code Total Time:      %8.4f seconds\n', time_ai);
fprintf('Speedup Factor:          %.2fx faster\n', time_manual / time_ai);


% -------------------------------------------------------------------
%  MANUAL SOLVER 
% -------------------------------------------------------------------
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

        [lower_discon, higher_discon] = discont_finder(2.3, 2.6, m);
        root_value = Bi_section_root_finder(lower_discon + 0.01,...
            higher_discon - 0.01, Bi, epsilon);
   % end
    
end

%to save run time, we made 2 txt files containing the discontinuites of 
%the function, which are the roots of the zeroeth bessel function, more 
%specifically the first 100 roots, and the first 10000, and to find the
%m-th root of our function, we accses the m-th and (m-1)-th 
%discontinuity and search for a root in those boundaries, assuring we 
%won't be looking around a discontinuity in the function, and won't
%need to check "m" roots of the function.
    % if m == 1
    %     root_value = Bi_section_root_finder(0.3, 2.4, Bi, epsilon);
    % elseif m <= 100
    %     J0_roots = readmatrix("roots_short.txt");
    %     %given the Biot numbers, a distance of 0.1 from the discontinuty
    %     % will suffice for every root. 
    %     higher_guess = J0_roots(m) - 0.1; 
    %     lower_guess = J0_roots(m-1) + 0.1;
    %     root_value = ...
    %         Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    % else
    %     J0_roots = readmatrix("roots_long.txt");
    %     %the readmatrix comand is slow for this size of a file, so we
    %     %preffered to split the roots files and open only the neccesary
    %     %one.
    % end
    % 
    % if m > 100 && m <= 10000
    %     higher_guess = J0_roots(m) - 0.1;
    %     lower_guess = J0_roots(m-1) + 0.1;
    %     root_value = ...
    %         Bi_section_root_finder(lower_guess, higher_guess, Bi, epsilon);
    % elseif m > 10000
    %     %if there is no other choise, we will continue mapping the
    %     %discontinuiues and find the root between the (m-1)-th and m-th
    %     %discontinuity points. this is obviuosly the longest and least
    %     %desired path to take. the "pi" distance is aprroximatly the
    %     %bessel function frequncy.

        % [lower_discon, higher_discon] = discont_finder(J0_roots(10000)...
        %     + 0.01, J0_roots(10000) + 0.01 + pi, m);

        function [func_value] = f(x, bi_num)
%The function will recive an x value, scalar or array, and a biot number,
%and will return the value of x multiplyed by J_1(x)/J_0(x) minus the biot
%number
    
    func_value = x.*(besselj(1,x)./besselj(0,x)) - bi_num;

end

function [x_prev, x_m] = discont_finder(xIn1, xIn2, m)
%This function will return the (m-1)-th and m-th roots of the
%zeroeth bessel function after recieving 2 initial values to find the
%first root between them, and the wanted m.
    
    %initilizing the parameters, and aux variables.
    flag = 0;
    sec_flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    x_prev = 0;
    x_m = 0;
    if m <= 10000
        m = m + 10000;
    end

    for i = 1:m-10000
        %the loop will run m-10000 times to return the m-th root while
        %constantly keeping the previous root.
        if i > (m - 10000 - 2)
                epsilon = 10^(-7);
                %improving the accuracy in the last roots, the ones that
                %will be used.
        else
            epsilon = 0.0001;
        end
        while(flag == 0)
            if counter > 2000
                error("Too many iterations, we have a problem")
                %preventing infinte loop, the counter is adding in the end
                %of it
            end
            
            if sec_flag == 0
                f_lower = besselj(0,lower);
                %saving the new evalution of the function each time the
                %lower boundry isn't changing.
            end
            sec_flag = 0;    
            x_new = (lower + upper)/2;
            f_xnew = besselj(0, x_new);
            %calculating the value of the new x, assesing its function
            %value.
            
            if abs(f_xnew) < epsilon
                if i > (m - 10000 - 2)
                    x_prev = x_m;
                    x_m = x_new;
                end
                flag = 1;
                %each root being identified, but only the relevant one
                %stored
                
            elseif f_lower*f_xnew < 0
                upper = x_new;
                sec_flag = 1;
                %if not close enough to a root, cutting by half the range
                %and checking for different signs of the function at the
                %boundries.
            elseif f_lower*f_xnew > 0
                lower = x_new;
%practically, the rest of those if statements are meaningless, it will
%never reach them.
            elseif f_lower == 0
                if i > (m - 10000 - 2)
                    x_prev = x_m;
                    x_m = lower;
                end
                flag = 1;
                %each root being identified, but only the relevant one
                %stored
            elseif f(upper, Bi) == 0
                if i > (m - 10000 - 2)
                    x_prev = x_m;
                    x_m = upper;
                end
                flag = 1;
                %each root being identified, but only the relevant one
                %stored
                
            end
            
            counter = counter + 1;
            
        end

        %initilizing again for the next root
        
        flag = 0;
        sec_flag = 0;
        counter = 0;
        
        lower = x_new + 0.1;
        upper = x_new + 0.1 + pi;
    end
end

function [root_value] = Bi_section_root_finder(xIn1, xIn2, Bi_num, epsilon)
%This function utelize the Bi section method for finding a root of a
%function in a limited range - the function will recive boundries, biot
%number and tolerance, and return an 'x' which its function value is within
%an epsilon from zero.

    %initilizing the parameters and aux variables
    flag = 0;
    sec_flag = 0;
    counter = 0;
    lower = xIn1;
    upper = xIn2;
    root_value = 0;

    if f(lower, Bi_num)*f(upper, Bi_num) > 0
        upper = lower + pi/2 + 0.1;
        %for cases the boundries weren't defined correctly, this will
        %provide a result, although less reliable one.
    end
            
    while(flag == 0)
        if counter > 10000
            error("Too many iterations, we have a problem")
            %prevents an endless loop
        end
        
        %to prevent it from running the besselj function each time, it will
        %do only if the lower x value had been changed.
        if sec_flag == 0
            f_lower = f(lower, Bi_num);
        end
        sec_flag = 0;
        
        x_new = (lower + upper)/2;
        %x_new will be the average of the prevoius 'x's.
        f_xnew = f(x_new, Bi_num);
        
        if abs(f_xnew) < epsilon
            %if the new x value isn't close enough, we will choose the x
            %that yield result with different sign than the current.
            root_value = x_new;
            flag = 1;
        elseif f_lower*f_xnew < 0
            upper = x_new;
            sec_flag = 1;
        else
            lower = x_new;
        end
        
        counter = counter + 1;
    end
end

