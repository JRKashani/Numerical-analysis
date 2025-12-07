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