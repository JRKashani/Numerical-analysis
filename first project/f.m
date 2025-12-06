function [func_value] = f(x, bi_num)

%    if isvector(x)
 %      func_value = zeros(1, length(x));
  %     for i = 1:length(x)
   %        func_value(i) = x(i)*(besselj(1,x(i))/besselj(0,x(i))) - bi_num;
    %   end
    %else
     %   func_value = x*(besselj(1,x)/besselj(0,x)) - bi_num;
    %end
    
    func_value = x.*(besselj(1,x)./besselj(0,x)) - bi_num;

end

