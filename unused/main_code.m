function [X]=main_code(Bi,m)
if isinteger(m) || ~isscalar(m)
    error("Irrelevant input.");
end
%if m>=10000
    upper=2.5;
    lower=0.1;
    epsilon=1e-2;
    flag=0;
    count=0;
    n=m;
    for i=(1:n)
        while flag==0
            count=count+1;        
            m_root=(upper+lower)/2;
            if count==10000
                error('too much eteration');
            elseif abs(besselj(0,m_root))<epsilon && besselj(0,lower)*besselj(0,m_root)<0
                flag=1;
            elseif besselj(0,lower)*besselj(0,m_root)<0
                upper=m_root;
            else
                lower=m_root;
            end
        end
        flag=0;
        count=0;
        upper=m_root+pi;
        lower=m_root;
    end
%end
X=cordway(Bi,lower,upper);
end
