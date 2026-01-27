function [outputArg1,outputArg2] = FE_solver(inputArg1,inputArg2)


    while
        u(i + 1, 1) = u(i, 1) + h;
        u(i + 1, 2) = u(i, 2) + h*f1(u(i, 1), u(i, 2), u(i, 3), u(i, 4));
        u(i + 1, 3) = u(i, 3) + h*f2(u(i, 1), u(i, 2), u(i, 3), u(i, 4));
        u(i + 1, 4) = u(i, 4) + h*f3(u(i, 1), u(i, 2), u(i, 3), u(i, 4));
    end
end