function [A, err_weights, b5, c] = RKF_Butcher_Tableau(~)
%The function would recieve no input arguments and will return a matrix and
%3 vectors of the parameters needed for the RKF method.

    %matrix A for 6 k's
    A = [0           0           0           0           0      0;
         1/4         0           0           0           0      0;
         3/32        9/32        0           0           0      0;
         1932/2197  -7200/2197   7296/2197   0           0      0;
         439/216    -8           3680/513   -845/4104    0      0;
        -8/27        2          -3544/2565   1859/4104  -11/40  0];

    %weight for 5th order solutions.
    b5 = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];

    %weight for 4th order solutions, for RKF it will be used only as error
    %control methos
    b4 = [25/216; 0; 1408/2565;  2197/4104;   -1/5;  0];

    %instead of computing the actual values of u5 and u4 and subtracting,
    %the weights are subtracted
    err_weights = b5 - b4;

    %Time nodes
    c  = [0; 1/4; 3/8; 12/13; 1; 1/2]; 
end