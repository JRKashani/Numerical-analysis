function zero_x = find_J0_zero_bisection(k, TOL, MAX_ITER)
% FINDS_J0_ZERO_BISECTION Finds the k-th positive zero of J0(x) using 
% the Bisection Method. This is used to define the starting bracket
% for the main root-finding problem.
%
% Inputs:
%   k - The zero number (k >= 1)
%   TOL - Tolerance for convergence
%   MAX_ITER - Maximum iterations
%
% Output:
%   zero_x - The k-th zero of J0(x).

    f_J0 = @(x) besselj(0, x);
    
    % --- Initial Bracket Selection and Justification ---
    % The k-th zero of J0(x), j_{0, k}, is approximately pi*(k - 1/4).
    % The distance between consecutive zeros is approximately pi.
    % We search for an interval [a, b] where J0(a) * J0(b) < 0.
    
    % Start the search near the (k-1)-th zero and increment until a sign change is found.
    a = pi * (k - 1); 
    b = pi * k;
    step = pi/2; % A safe search step size

    % Refine the bracket until f_J0(a) and f_J0(b) have opposite signs
    while f_J0(a) * f_J0(b) >= 0
        a = b;
        b = b + step;
        % Safety limit for search (ensures termination for very large k)
        if b > pi * k + 5*pi 
            error('Could not bracket the k-th zero of J0(x). Check the initial approximation logic.');
        end
    end
    
    % --- Bisection Method (Guaranteed Convergence) ---
    for iter = 1:MAX_ITER
        mid = (a + b) / 2;
        f_mid = f_J0(mid);
        
        if abs(f_mid) < TOL % Convergence criterion
            zero_x = mid;
            return;
        end
        
        if f_J0(a) * f_mid < 0 % Root is in [a, mid]
            b = mid;
        else                   % Root is in [mid, b]
            a = mid;
        end
    end
    
    warning('Bisection for J0 zero did not converge.');
    zero_x = (a + b) / 2; % Return the best approximation
end