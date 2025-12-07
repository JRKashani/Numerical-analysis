function [f, df] = bessel_func_product(x, Bi)
    % Calculate Bessel functions of first kind, orders 0 and 1
    J0 = besselj(0, x);
    J1 = besselj(1, x);
    
    % The function value: x*J1 - Bi*J0
    f = x .* J1 - Bi .* J0;
    
    % The derivative: x*J0 + Bi*J1
    df = x .* J0 + Bi .* J1;
end
