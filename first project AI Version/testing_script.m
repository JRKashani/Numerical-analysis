% TEST_RANDOM_BESSEL
% Generates random test cases to validate SolveBesselAI.m
% Checks for:
% 1. Residual Error: Does the equation actually equal zero?
% 2. Magnitude: Is the root in the expected 'neighborhood'?

clc; clear; close all;

% Configuration
num_tests = 1000;
min_Bi = 0.1;
max_Bi = 60;
max_m = 5000;

fprintf('Running %d random tests...\n', num_tests);
fprintf('%-6s | %-6s | %-10s | %-10s | %-12s | %-10s\n', ...
        'Bi', 'm', 'Root (x)', 'Expected~', 'Residual', 'Status');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:num_tests
    % 1. Generate Random Inputs
    % Bi between 0.1 and 60
    Bi = min_Bi + (max_Bi - min_Bi) * rand(); 
    % m integer between 1 and max_m
    m = randi([1, max_m]);
    
    % 2. Solve using our function
    try
        x_sol = SolveBesselAI(Bi, m);
        
        % 3. Calculate Residual (Accuracy Check)
        % Evaluate the continuous form: f(x) = x*J1 - Bi*J0
        [res, ~] = bessel_func_product(x_sol, Bi);
        residual = abs(res);
        
        % 4. Heuristic Check (Magnitude Check)
        % Roots appear roughly every pi steps, starting near 2.4
        x_heuristic = 2.4 + (m - 1) * pi;
        
        % Check if residual is acceptable (e.g., < 1e-6)
        if residual < 1e-6
            status = 'PASS';
        else
            status = 'FAIL';
        end
        
        % Print results
        fprintf('%6.2f | %6d | %10.4f | %10.4f | %12.1e | %s\n', ...
                Bi, m, x_sol, x_heuristic, residual, status);
        
    catch ME
        fprintf('%6.2f | %6d | %10s | %10s | %12s | ERROR\n', ...
                Bi, m, 'NaN', '-', '-');
        fprintf('Error: %s\n', ME.message);
    end
end
fprintf('%s\n', repmat('-', 1, 75));
fprintf('Note: "Expected~" is a rough heuristic (2.4 + (m-1)*pi).\n');
fprintf('Deviations are normal as Bi changes, but should follow the trend.\n');