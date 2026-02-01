function [t, U] = RK4_solver(func, t_span, u0, para)
%The function will recieve a parameter struct, a function and initial
%conditions, and will return the altered struct and a matrix of values - a
%row vector for each variable.
   
    h = para.delta; %unpacking the step size
    if isscalar(t_span)
        t_span = [0, t_span];
    end
    t = t_span(1) : h : t_span(end); %Generate time vector

    N = length(t);

    if length(u0) ~= para.ode_rank
        error("Initial values are not in the required order and amount");
    end

    U = zeros(para.ode_rank, N); %each stat get a row, each data 
    % point (time-value) gets a coloumn.
    U(:, 1) = u0(:); %

    for i = 1:N-1
        uk = U(:, i);
        k1 = func(t(i),      uk); %dervs at the original point, u_k
        k2 = func(t(i)+ h/2, uk + ((h/2) * k1)); %dervs at the point 
                %which will get to after a half step in f1 direction
        k3 = func(t(i)+ h/2, uk + ((h/2) * k2)); %dervs at the point 
                %which will get to after a half step in f2 direction
        k4 = func(t(i)+ h,   uk +  (h    * k3)); %dervs at the point 
                %which will get to after a full step in f3 direction
        
        U(:, i + 1) = uk + ((h/6) * (k1 + 2*k2 + 2*k3 + k4)); %the new 
        % stats-vector is the original point plus a weighted average of
        % the dervs multiplyed by the time step        
    end
end