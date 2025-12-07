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