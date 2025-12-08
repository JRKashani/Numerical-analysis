% --- Helper: Robust Bisection ---
function [root, iter] = run_bisection(func, a, b, tol, max_iter)
    [fa, ~] = func(a);
    [fb, ~] = func(b);
    
    % Final safety check before burning iterations
    if fa * fb > 0
         error('Bisection failed: Bracket invalid (same signs)');
    end
    
    for iter = 1:max_iter
        c = (a + b) / 2;
        [fc, ~] = func(c);
        
        if abs(fc) < tol || (b - a)/2 < tol
            root = c;
            return;
        end
        
        if fc * fa > 0
            a = c; fa = fc;
        else
            b = c; % fb = fc;
        end
    end
    root = (a + b) / 2;
end