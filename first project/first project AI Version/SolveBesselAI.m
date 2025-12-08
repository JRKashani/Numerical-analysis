function root = SolveBesselAI(Bi, m)
% SOLVEBESSELAI Finds the m-th positive root of x*J1(x) = Bi*J0(x).
%
% Inputs:
%   Bi : Biot number (0.1 <= Bi <= 60)
%   m  : Index of the root (1, 2, 3...)
%
% Strategy:
%   1. Scan for the bracket [x_min, x_max] containing the m-th root.
%   2. Initial Guess: Secant method approximation (smarter than midpoint).
%   3. Refine: Constrained Newton-Raphson (must stay in bracket).
%   4. Fallback: Robust Bisection if Newton fails or leaves bracket.

    % --- Configuration ---
    tol = 1e-8;      % Convergence tolerance
    max_iter = 100;  % Iteration limit
    
    % --- Step 1: Input Validation ---
    if Bi < 0.1 || Bi > 60
        error('Input Error: Bi must be between 0.1 and 60.');
    end
    if m <= 0 || mod(m, 1) ~= 0
        error('Input Error: m must be a positive integer.');
    end

    % --- Step 2: Define Function Handle ---
    target_func = @(x) bessel_func_product(x, Bi);

    % --- Step 3: Robust Scanning ---
    % Step size: pi/20 gives ~15-20 points per oscillation period (pi),
    % ensuring we don't skip roots even for high m.
    x_curr = 1e-4; 
    step_size = pi / 20; 
    roots_found = 0;
    
    % Limit: Roots are roughly at (m + 0.75)*pi. We add buffer.
    search_limit = (m + 1) * pi + 10; 
    
    x_min = NaN;
    x_max = NaN;
    [f_next, ~] = target_func(x_curr);
    bracket_found = false;

    while ~bracket_found
        x_next = x_curr + step_size;
        
        f_curr = f_next;
        [f_next, ~] = target_func(x_next);
        
        % Check for sign change
        if f_curr * f_next < 0
            roots_found = roots_found + 1;
            if roots_found == m
                x_min = x_curr;
                x_max = x_next;
                bracket_found = true;
            end
        end
        
        x_curr = x_next;
        
        if x_curr > search_limit
            error('Scan failed: Root #%d not found within x < %.1f', m, search_limit);
        end
    end

    % --- Step 4: Hybrid Solver ---
    try
        % Smart Initial Guess (Secant Approximation)
        % This places x0 closer to the root than a simple midpoint.
        [fa, ~] = target_func(x_min);
        [fb, ~] = target_func(x_max);
        
        % Safety: Ensure bracket is valid before starting
        if fa * fb > 0
            error('Invalid bracket: endpoints have same sign');
        end
        
        % Linear interpolation for x0
        x0 = x_min - fa * (x_max - x_min) / (fb - fa);

        % Attempt Constrained Newton-Raphson
        [root, ~] = run_newton(target_func, x0, tol, max_iter,...
            x_min, x_max);
        
    catch
        % Fallback to Bisection
        % (Triggered if Newton diverges, leaves bracket, or hits
        % zero derivative)
        [root, ~] = run_bisection(target_func, x_min, x_max, tol,...
            max_iter*10);
    end
end
