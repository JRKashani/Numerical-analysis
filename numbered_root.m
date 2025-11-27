function [root_value] = numbered_root(Bi, m)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if ~isinteger(m) || ~isscalar(m)
        error("Irrelevant input.");
    elseif m < 10000
        vector_of_J0_roots = readmatrix("ROOTS.txt", 'NumVariables', m);
    else
        vector_of_J0_roots = readmatrix("ROOTS.txt");
    end

    
    


end