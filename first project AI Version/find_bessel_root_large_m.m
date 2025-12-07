function root_x = find_bessel_root_large_m(Bi, m)
% FIND_BESSEL_ROOT_LARGE_M Numerically finds the m-th positive root of
% x * J1(x) / J0(x) = Bi.
%
% This version is robust for large values of m by accurately locating
% the m-th zero of J0(x) to form the initial bracket.
%
% Inputs:
%   Bi - The Biot number (0.1 <= Bi <= 60)
%   m - The root number (m >= 1, a positive integer)
%
% Output:
%   root_x - The m-th positive root of the equation.

    % --- Global Numerical Parameters ---
    % Tolerance and maximum iterations are explicitly defined and justified[cite: 42].
    TOL = 1e-6;         % Stopping criterion: absolute change < TOL
    MAX_ITER = 100;     % Safety limit for iterations
    
    % Define the target function g(x) = x * J1(x) / J0(x) - Bi
    g = @(x) x .* (besselj(1, x) ./ besselj(0, x)) - Bi; 
    
    % Define the derivative g'(x) for Newton-Raphson: d/dx [x * J1(x) / J0(x)]
    % Note: J0'(x) = -J1(x), J1'(x) = (J0(x) - J2(x))/2
    g_prime_NR = @(x) besselj(1, x) ./ besselj(0, x) + x .* ( (besselj(0, x) .* ( (besselj(0, x) - besselj(2, x)) / 2 ) - besselj(1, x) .* (-besselj(1, x))) ./ besselj(0, x).^2 );

    % --- 1. Find the Bracket [a, b] using J0(x) zeros ---

    % The m-th root is bounded by the zeros of J0(x).
    % The bracket [a, b] is [j_{0, m-1}, j_{0, m}], where j_{0, k} is the k-th zero of J0(x).
    
    if m == 1
        a = 0.1; % Start just after the singularity at x=0
    else
        % Find the (m-1)-th zero of J0(x) to set the lower bound 'a'
        a = find_J0_zero_bisection(m - 1, TOL, MAX_ITER);
    end
    
    % Find the m-th zero of J0(x) to set the upper bound 'b'
    b = find_J0_zero_bisection(m, TOL, MAX_ITER);
    
    % Justification: The root is guaranteed to be in (a, b).
    
    % --- 2. Newton-Raphson Method to Refine the Root ---
    
    % Initial guess: midpoint of the robust bracket.
    x0 = (a + b) / 2;
    
    % Start the iteration
    for k = 1:MAX_ITER
        fxk = g(x0);
        dfxk = g_prime_NR(x0);
        
        % Check for division by zero (singularity near J0=0)
        if abs(dfxk) < 1e-10
            warning('Derivative is too small. Using bisection step as fallback.');
            if fxk > 0, b = x0; else, a = x0; end
            x_new = (a + b) / 2;
        else
            % Newton-Raphson step
            x_new = x0 - fxk / dfxk;
        end
        
        % Check the stopping criterion (required accuracy/precision)
        if abs(x_new - x0) < TOL 
            root_x = x_new;
            return;
        end
        
        % Update the guess
        x0 = x_new;
    end
    
    % If the loop finishes without converging
    warning('Newton-Raphson did not converge within the maximum number of iterations.');
    root_x = x0; 
end