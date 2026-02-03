function [df_dist_from_one] = Blasius_Wrapping(ddf_at_zero)
% The function used to solve the Blasius boundary equation with the given
% initial values. the only input is f''(0), and the output is the limit in
% inf for f'(eta)

    %% Parameters
    para.ode_rank  = 3;
    %para.delta     = 0.1;     % step size for RK4, legacy
    para.epsilon   = 10^(-8);   % tolerance for RKF
    para.h_init    = 0.005;     % very small, but will be adapted swiftly
    para.h_max     = 0.1;      % to prevent the solver miss changes.
    para.max_steps = 20000;     %machine dependent, theoretically
    para.safety_factor = 0.841; %0.5^0.25, will elaborate at the notebook
    

    t_span = [0, 10]; %the function won't change much after 10, but the
    % cost of the safety factor was negligable.
    u0 = [0; 0; ddf_at_zero]; %initial values vector : f(0) = 0; f'(0) = 0;
    % f''(0) = given value.
    
    func = @(t, u) [ u(2); u(3); (-0.5) * u(1) * u(3)]; %The blasius 
    % boundary equation, linearised.

    [eta, U] = RKF_solver(func, t_span, u0, para); %the eta values are 
    % omitted, given the arent relevant to the assignment
   
    df_dist_from_one = U(2, end) - 1; %regardless to noise, the last
    % computed value of f'(eta) (highest eta) is the closest to the limit
end


% function [df_dist_from_one, dg] = Blasius_Wrapping(ddf_at_zero)
% % The function used to solve the Blasius boundary equation simoltanusly to
% % the function that evaluate d(f'(eta))/d(f''(0)) in each data point,
% % serving as the needed derivitive for the Newton-Raphson root finding
% % method. The function accept the initial value f''(0) and return the
% % distance of f'(inf) from 1, and the gradient of said distance inrespect
% % to the input initial value
% 
%     %% Parameters
%     para.ode_rank  = 6; %3 variables 
%     %para.delta     = 0.1;     % step size for RK4
%     para.epsilon   = 10^(-8);    % tolerance for RKF
%     para.h_init    = 0.005;
%     para.h_max     = 0.01;
%     para.max_steps = 20000;
%     para.safety_factor = 0.841;
%     %delta_t = sqrt(para.epsilon);
% 
%     t_span = [0, 12];
%     u0 = [0; 0; ddf_at_zero; 0; 0; 1];
%     %u1 = [0; 0; ddf_at_zero + delta_t; 0; 0; 1];
% 
%     func = @(t, u) [ u(2); u(3); (-0.5) * u(1) * u(3);...
%         u(5); u(6); (-0.5) * (u(4)*u(3) + u(6)*u(1)) ];
%     [~, U] = RKF_solver(func, t_span, u0, para);
%     %[~, U2] = RKF_solver(func, t_span, u1, para);
% 
%     df_dist_from_one = U(2, end) - 1;
%     %dg = (U2(2, end) - U(2, end))/delta_t;
%     dg = U(5, end); % + 2 * U(6, end);
%     %ddg = U(6, end);
% end
% 

