function root = SolveBesselAI(Bi, m)
% SOLVEBESSELAI Main solver for the Bessel project (AI-assisted).
%   root = SolveBesselAI(Bi, m) finds the m-th positive root of the
%   equation x*J1(x) = Bi*J0(x).
%
%   Strategy:
%   1. Scan for the interval containing the m-th root.
%   2. Attempt Newton-Raphson (Fastest, O(h^2)).
%   3. Fallback to Bisection if Newton fails (Robust, O(1/2^n)).

    % --- Configuration ---
    tolerance = 1e-8;      % High precision
    max_iter = 100;        % Sufficient for Newton
    
    % --- Step 1: Input Validation ---
    validate_inputs(Bi, m);

    % --- Step 2: Define the Function Handle ---
    % We use the continuous form: f(x) = x*J1(x) - Bi*J0(x)
    % This handle encapsulates 'Bi' so inner functions don't need it.
    target_func = @(x) bessel_func_product(x, Bi);

    % --- Step 3: Bracket the Root ---
    % Find the interval [x_min, x_max] where the root exists
    [x_min, x_max] = scan_for_bracket(target_func, m);

    % --- Step 4: Hybrid Numerical Solution ---
    try
        % ATTEMPT 1: Newton-Raphson
        % We choose the midpoint as the best unbiased guess
        initial_guess = (x_min + x_max) / 2;
        
        [root, ~] = solve_newton(target_func, initial_guess, ...
                                 tolerance, max_iter);
                                 
    catch
        % ATTEMPT 2: Bisection (Fallback)
        % Only runs if Newton fails (e.g., derivative is zero)
        % warning('Newton failed. Switching to Bisection.'); % Optional
        
        [root, ~] = solve_bisection(target_func, x_min, x_max, ...
                                    tolerance, max_iter * 10); 
                                    % Bisection needs more iterations
    end
end

function [root, iter] = solve_newton(func_handle, x0, tol, max_iter)
% Inputs:
%   func_handle: Handle to the function returning [f, df]
%   x0:          Initial guess
%   tol:         Tolerance for stop condition
%   max_iter:    Maximum allowed iterations

    x_current = x0;

    for iter = 1:max_iter
        [f, df] = func_handle(x_current);

        % Stop Condition 1: Function value is close enough to zero
        if abs(f) < tol
            root = x_current;
            return;
        end

        % Safety: Avoid division by zero if derivative vanishes
        if abs(df) < 1e-10 
            error('Derivative is too close to zero. Newton method failed.');
        end

        % Newton Step
        x_next = x_current - (f / df);

        % Stop Condition 2: Step size is smaller than tolerance (convergence)
        if abs(x_next - x_current) < tol
            root = x_next;
            return;
        end

        x_current = x_next;
    end

    error('Newton-Raphson method did not converge within max_iter.');
end

function [root_approximation, iterations_used] = solve_bisection(...
    target_function, lower_bound, upper_bound, tolerance, max_iterations)
% Inputs:
%   target_function: Handle to the function returning [value, derivative]
%   lower_bound:     The start of the search interval
%   upper_bound:     The end of the search interval
%   tolerance:       The convergence threshold
%   max_iterations:  Limit to prevent infinite loops

    % Evaluate the function at the boundaries
    [val_at_lower, ~] = target_function(lower_bound);
    [val_at_upper, ~] = target_function(upper_bound);

    % Validate that the root is bracketed (signs must be opposite)
    if sign(val_at_lower) == sign(val_at_upper)
        error(['Error: The interval [lower_bound, upper_bound] does' ...
            ' not bracket a root (signs are identical).']);
    end

    for current_iter = 1:max_iterations
        % Calculate the midpoint
        midpoint = (lower_bound + upper_bound) / 2;
        [val_at_mid, ~] = target_function(midpoint);

        % Stop Condition 1: Function value is sufficiently close to zero
        % Stop Condition 2: The interval has shrunk below the tolerance
        current_interval_width = (upper_bound - lower_bound);
        if abs(val_at_mid) < tolerance || (current_interval_width / 2)...
                < tolerance
            root_approximation = midpoint;
            iterations_used = current_iter;
            return;
        end

        % Narrow the interval based on the sign change
        if sign(val_at_mid) == sign(val_at_lower)
            % The root is in the right half
            lower_bound = midpoint;
            val_at_lower = val_at_mid;
        else
            % The root is in the left half
            upper_bound = midpoint;
            % val_at_upper = val_at_mid; % Optional update
        end
    end
    
    error(['Bisection method failed to converge within the maximum' ...
        ' number of iterations.']);
end

function validate_inputs(Biot_number, root_index_m)
% Inputs:
%   Biot_number:  Dimensionless number representing heat transfer ratio
%   root_index_m: The index of the positive root to find (1, 2, 3...)

    % 1. Validate Biot Number Range
    % Source constraint: 0.1 <= Bi <= 60 [cite: 252]
    if Biot_number < 0.1 || Biot_number > 60
        error(['Input Error: Biot number (Bi) must be between 0.1 ' ...
            'and 60. You provided: %.4f'], Biot_number);
    end

    % 2. Validate Root Index (m)
    % Source constraint: m is a natural number [cite: 255]
    % Check if m is positive and if it is an integer
    if root_index_m <= 0 || mod(root_index_m, 1) ~= 0
        error(['Input Error: The root index (m) must be a natural' ...
            ' number (positive integer). You provided: %f'], root_index_m);
    end
end

function [x_min, x_max] = scan_for_bracket(func_handle, target_m)
% SCAN_FOR_BRACKET Scans for the m-th root interval.
% Inputs:
%   func_handle: The continuous function (x*J1 - Bi*J0)
%   target_m:    The index of the root we want (1st, 2nd, etc.)
% Output:
%   [x_min, x_max]: The interval containing the m-th root.

    % Start slightly above 0 to avoid trivial solutions/singularities
    x_current = 1e-4; 
    step_size = 0.1;
    
    roots_found_count = 0;
    bracket_found = false; % Flag to control the loop
    
    % Initialize outputs to avoid errors if loop fails
    x_min = NaN;
    x_max = NaN;
    
    while ~bracket_found
        x_next = x_current + step_size;
        
        [f_curr, ~] = func_handle(x_current);
        [f_next, ~] = func_handle(x_next);
        
        % Check for sign change
        if sign(f_curr) ~= sign(f_next)
            roots_found_count = roots_found_count + 1;
            
            % Check if this is the target m-th root
            if roots_found_count == target_m
                x_min = x_current;
                x_max = x_next;
                bracket_found = true; % Set flag to exit loop
            end
        end
        
        % Move to next step
        x_current = x_next;
        
        % Safety exit: Prevent infinite loops if inputs are bad
        if x_current > 1000 
            error('Exceeded search limit. Root m=%d not found.', target_m);
        end
    end
end
