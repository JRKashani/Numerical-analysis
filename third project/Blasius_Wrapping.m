function [df_dist_from_one, dg] = Blasius_Wrapping(ddf_at_zero)

    
    %% Parameters
    para.ode_rank  = 6;
    %para.delta     = 0.1;     % step size for RK4
    para.epsilon   = 10^(-8);    % tolerance for RKF
    para.h_init    = 0.05;
    para.h_max     = 0.1;
    para.max_steps = 20000;
    para.safety_factor = 0.841;
    %delta_t = sqrt(para.epsilon);

    t_span = [0, 12];
    u0 = [0; 0; ddf_at_zero; 0; 0; 1];
    %u1 = [0; 0; ddf_at_zero + delta_t; 0; 0; 1];
    
    func = @(t, u) [ u(2); u(3); (-0.5) * u(1) * u(3);...
        u(5); u(6); (-0.5) * (u(4)*u(3) + u(6)*u(1)) ];
    [~, U] = RKF_solver(func, t_span, u0, para);
    %[~, U2] = RKF_solver(func, t_span, u1, para);
   
    df_dist_from_one = U(2, end) - 1;
    %dg = (U2(2, end) - U(2, end))/delta_t;
    dg = U(5, end); % + 2 * U(6, end);
    %ddg = U(6, end);
end
