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