function [root_value] = numbered_root_3ver(Bi, m)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if isinteger(m) || ~isscalar(m)
        error("Irrelevant input.");
    elseif m < 10000
        load('m_valious.mat','m_mat');
    end
    if m==1
    upper_boundry = m_mat(1,m) - 0.1;
    x=upper_boundry;
    else
    lower_boundry = m_mat(1,m-1)-0.1;
    x=lower_boundry;
    end

    counter = 0;
    func_value=1;
    epsilon = 1e-5; % Define a small tolerance for convergence
    
    while true
    if abs(func_value) < epsilon
      root_value = x;
      break;
    elseif counter > 10000
       error(too_many_loops_text);
    end
    %in case of an unresoled question, the counter will prevent from
    %infinte loop
    counter = counter + 1;
    %getting the bessel functions values into compact form, for improved
    %readiness
    J0 = besselj(0, x);
    J1 = besselj(1, x);
    J2 = besselj(2, x);
    %calculationg the function and derived values:
    func_value = x*(J1/J0) - Bi;
    deriv_value = (J1/J0)*(1+(x*(J1/J0))) - x/2 - J2/(2*J0);
    %the Newton-Raphson method core
    y = x - (func_value/deriv_value);
    x = y; 
    end
end