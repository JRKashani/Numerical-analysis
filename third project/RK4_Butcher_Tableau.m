function [A, b4, c] = RK4_Butcher_Tableau(~)
%The function would recieve no input arguments and will return a matrix and
%2 vectors of the parameters needed for the RK4 method

    A = [ 0     0     0     0;
          1/2   0     0     0;
          0     1/2   0     0;
          0     0     1     0 ];

    b4 = [1/6; 1/3; 1/3; 1/6];

    c  = [0; 1/2; 1/2; 1];
end