%this isn't a function, but a script meant to be run once. we are adding
%this to the work to clearify how we worked, this file won't be use by any
%other function in the program.

%initilizers - the numbers are derived from plotting the function and
%finding the roots and the discontinuitues.
lower = 0.1;
upper = pi + 0.2;
counter = 0;
flag = 0;
sec_flag = 0;
epsilon = 1e-8;
%naturally, there is no need for such accuracy, but given it will run once
%or twice, there is no need for time efficency.

%allocating the memory for the arrays so they won't change size each
%itiration
J0_roots_short = zeros(1, 100);
J0_roots_long = zeros(1, 10000);

for i = 1:100
    while(flag == 0)
        if counter > 2000
            error("Too many iterations, we have a problem")
        end
        
        %to prevent it from running the besselj function each time, it will
        %do only if the lower x value had been changed.
        if sec_flag == 0
            f_lower = besselj(0,lower);
        end
        sec_flag = 0;

        x_new = (lower + upper)/2;
        %x_new will be the average of the prevoius 'x's.
        f_xnew = besselj(0, x_new);

        if abs(f_xnew) < epsilon
            %success testing, if pass, will raise a flage and the while
            %loop will end.
            J0_roots_short(i) = x_new;
            flag = 1;
        elseif f_lower*f_xnew < 0
            %if the new x value isn't close enough, we will choose the x
            %that yield result with different sign than the current. 
            upper = x_new;
            sec_flag = 1;
        else
            lower = x_new;
        end

        counter = counter + 1;
        %the counter prevent endless loops.

    end

    flag = 0;
    counter = 0;
    lower = x_new + 0.1;
    upper = x_new + 0.1 + pi;
    %when moving to the next root, we initilise the flag, counter and
    %assigning new boundries, derived from the pi frequncy of the function.
end

writematrix(J0_roots_short, 'roots_short.txt');
%storing as txt file

J0_roots_long(1:100) = J0_roots_short;
%saving a short second

lower = J0_roots_long(100) + 0.1;
upper = J0_roots_long(100) + pi + 0.1;
%initlizing the boundries again.

for i = 101:10000
    %very similar loop to the previous one, so no need for commenting.
    while(flag == 0)
        if counter > 2000
            error("Too many iterations, we have a problem")
        end
        
        f_lower = besselj(0,lower);
        
        x_new = (lower + upper)/2;
        f_xnew = besselj(0, x_new);
        
        if abs(f_xnew) < epsilon
            J0_roots_long(i) = x_new;
            flag = 1;
        elseif f_lower*f_xnew < 0
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