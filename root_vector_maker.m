
lower = 0.1;
upper = 2.5;
J0_roots = zeros(1,10000);

for i = (1:10000)
    J0_roots(i) = root_finder_J0(lower,upper);
    lower = J0_roots(i) + 0.1;
    upper = J0_roots(i) + pi+1;
    if J0(lower)*J0(upper) > 0
        disp(lower)
        disp(upper)
        error("maybe it's not pi");
    end
end