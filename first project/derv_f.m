function [derv_value] = derv_f(xIn)
    
    a = besselj(1, xIn)/besselj(0, xIn);
    b = xIn * ((besselj(1, xIn)/besselj(0, xIn))^2);
    c = besselj(2, xIn)/(2*besselj(0, xIn));
    
    derv_value = a + b - c - (xIn/2);
end

