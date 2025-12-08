function answer=cheking_answer(m)
g=0;
flag =0;

for i=1:m
    while flag==0
        g=g+0.1;
        x=project1(g,i);
        b=x*(besselj(1,x)/besselj(0,x));
        if g~=b
            disp(b);
            disp(g);
            disp(i);
            error('there is a problem in hear!');
        end
         if g==60
            flag=1;
        end
    end
    flag=0;
    g=0;
end
answer=1;
end
