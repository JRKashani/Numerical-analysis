function derv_f = simp_Blasius_dervs(x, func_and_dervs)
%The function input is eta and the corresponding state-vector, in this case
%3-element long vector, and return the derivitives of the 3 functions at
%the point eta
    derv_f = [
        func_and_dervs(2);
        func_and_dervs(3);
        func_and_dervs(1) * func_and_dervs(3) * (-0.5)
              ];
end