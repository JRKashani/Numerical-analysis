function [f, df] = bessel_func_product(x, Bi)
% BESSEL_FUNC_PRODUCT Evaluates the continuous Bessel function and 
% its derivative.
%
% Equation: f(x) = x*J1(x) - Bi*J0(x)
%
% Derivative Derivation (Product Rule + Recurrence Relations):
%   Term 1: d/dx[x*J1] = J1 + x*J1'
%   Term 2: d/dx[-Bi*J0] = -Bi*J0' = Bi*J1
%   Identity used: J1' = (J0 - J2) / 2  [Source: Project PDF]
%   Result: f'(x) = (1 + Bi)*J1 + (x/2)*(J0 - J2)

    % Calculate Bessel functions of first kind, orders 0, 1, 2
    J0 = besselj(0, x);
    J1 = besselj(1, x);
    J2 = besselj(2, x);
    
    % The function value: x*J1 - Bi*J0
    f = x .* J1 - Bi .* J0;
    
    % The corrected derivative based on your derivation
    df = (1 + Bi) .* J1 + (x ./ 2) .* (J0 - J2);
end
