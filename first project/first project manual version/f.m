function [func_value] = f(x, bi_num)
%The function will recive an x value, scalar or array, and a biot number,
%and will return the value of x multiplyed by J_1(x)/J_0(x) minus the biot
%number
    
    func_value = x.*(besselj(1,x)./besselj(0,x)) - bi_num;

end

