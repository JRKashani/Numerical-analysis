function [para, t, U] = RK4_solver(func, t_span, u0, para)
%The function will recieve a parameter struct, a function and initial
%conditions, and will return the altered struct and a matrix of values - a
%row vector for each variable.
   
    h = para.delta;
    t = [0 : h : t_span];

    U = zeros(length(t), para.ode_rank);
    uk = u0;
    U(1, :) = u0;
    
    N = length(t);

    for i = 1:N-1
        f1 = func(t(i),      uk); %dervs at the original point, u_k
        f2 = func(t(i)+ h/2, uk + ((h/2) * f1)); %dervs at the point 
                %which will get to after a half step in f1 direction
        f3 = func(t(i)+ h/2, uk + ((h/2) * f2)); %dervs at the point 
                %which will get to after a half step in f2 direction
        f4 = func(t(i+1),    uk +  (h    * f3)); %dervs at the point 
                %which will get to after a full step in f3 direction

        U(i + 1, :) = uk + ((h/6) * (f1 + 2*f2 + 2*f3 + f4)); %the new 
        % stats-vector is the original point plus a weighted average of
        % the dervs multiplyed by the time step
        uk = U(i + 1, :); %appending the current state to the stats matrix
    end
    

end