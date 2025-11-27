function [func_value] = f(x, Bi)

    func_value = x*(besselj(1,x)/besselj(0,x)) - Bi;

end