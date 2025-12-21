%Jacobi matrix solution.
tic;

n = length(b);

x = zeros(n,1);
x_new = x;
epsilon = 1e-3;
flag = 0;
counter = 0;

L = tril(A,-1);
D = diag(A);
U = triu(A,1);
B=-1./D.*(L+U);
c=(1./D).*b;

while flag == 0
    counter = counter + 1;
    if counter > 1000
        error('Too many iterations');
    end

    x_new = B*x+c;
    
    if max(abs(x_new-x)./abs(x_new)) < epsilon
        flag = 1;
    end
    
    x = x_new;
end

toc;

%GS non-matrix solution.
tic

n = length(b); %Number of equations

x = zeros(1,n); %Initial guess.
x_new = x; %Next step
epsilon = 1e-3;
flag = 0;
counter = 0;
while flag == 0
    counter = counter + 1;
    if counter > 1000
        error('Too many iterations');
    end
    
    %GS Step
    for i = 1:n
        x_new(i) = (-sum(A(i,1:i-1).*x_new(1:i-1)) - ...
            sum(A(i,i+1:end).*x(i+1:end)) + b(i))/A(i,i);
    end
    
    %Stop condition
    if max(abs(x_new-x)./abs(x_new)) < epsilon
        flag = 1;
    end
    
    x = x_new;
end

toc