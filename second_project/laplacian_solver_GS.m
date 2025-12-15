% function [temp_mat] = laplacian_solver_GS(p)
% 
%     flag = 0;
%     counter = 0;
%     x_old = p.temp_mat(:);
%     x_new = x_old;
%     b = p.boundaries(:);
%     A = p.laplacian_mat;
%     n_rows = p.mat_size_1;
%     boundary = p.the_boundary;
% 
%     while flag == false
%         if counter > 1000
%             error("Too many iterations, something went wrong!");
%         end
%         counter = counter + 1;
% 
%         % x_new(1) =   b(1)/A(1, 1) - sum(A(1, 2 : end)' .* x_old(2 : end));
%         % x_new(end) = b(end)/A(end, end) - ...
%         %     sum(A(end, 1:(end-1))' .* x_new(1:(end-1)));
%         for i = n_rows + 1 : p.n_parameters - n_rows
%             %future_x = A(i, 1 : (i-1))'    .* x_new(1 : (i-1));
%             %past_x   = A(i, (i+1) : end)'  .* x_old((i+1) : end);
%             %x_new(i) = (b(i) - sum(future_x) - sum(past_x)) / A(i,i);
%             if ~boundary(i)
%               x_new(i) = (x_new(i - 1) + x_new(i - n_rows) +...
%                           x_old(i + 1) + x_old(i + n_rows)) / 4;
%             end
%         end
% 
%         if max(abs(x_new - x_old)) < p.epsilon
%             flag = true;
%         end
%         x_old = x_new;
%     end
%     temp_mat = reshape(x_new, p.mat_size_1, p.mat_size_2);
% end

function [p] = laplacian_solver_GS(p)

    temp_mat = p.temp_mat;
    lap_mat = p.laplacian_mat;
    N = p.n_parameters;
    flag = false;
    counter = 0;

    L = tril(lap_mat, -1);
    D = spdiags(diag(lap_mat), 0, N, N);
    U = triu(lap_mat, 1);
    LD = L+D;

    x = temp_mat(:);
    x_prev = zeros(length(x), 1);

    b = zeros(N, 1);
    b(p.the_boundary) = temp_mat(p.the_boundary);

    while flag == false
        
        c = b - U * x;
        for i = 1:N
            x(i) = (c(i) - (LD(i, 1:(i-1)) * x(1:i-1))) / LD(i , i);
        end

        counter = counter + 1;
        if counter > 10000
            error("Too much!");
        end

        if max(abs(x_prev - x)) < p.epsilon
            flag = true;
        end

        x_prev = x;
    end
    p.temp_mat = reshape(x, p.mat_size_1, p.mat_size_2);
end