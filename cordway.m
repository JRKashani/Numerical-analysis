function Ans=cordway(Bi,x0,x1)
epsilon=1e-4;
x_pre=x0;%x_n-1
x_corent=x1;%x_n
triger=0;
count=0;
while triger==0
    if count>1000000
        error('Too many iterations');
    else
        x_next=(f(x_corent,Bi)*x_pre-f(x_pre,Bi)*x_corent)/(f(x_corent,Bi)-f(x_pre,Bi));
    end
        if abs((g(x_next))-Bi)<epsilon
            triger=1;
        else
            count=count+1;
            x_pre=x_corent;
            x_corent=x_next;
        end
end
Ans=x_next;
end