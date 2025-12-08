function [nat_or_not] = is_natural(x)
    
    if isnumeric(x) && x > 0 && x == floor(x)
        nat_or_not = true;
    else
        nat_or_not = false;
end