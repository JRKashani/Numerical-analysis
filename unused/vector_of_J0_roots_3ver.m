
upper = 2.5;
J0_roots = zeros(1,100);

for i = (1:100)
    J0_roots(i) = root_finder_J0(lower,upper);
    lower = J0_roots(i) + 0.1;
    upper = J0_roots(i) + pi+1;
    if J0(lower)*J0(upper) > 0
        disp(lower)
        disp(upper)
        error("maybe it's not pi");
    end
    writematrix(J0_roots,'M_tab.txt');
   
end
 filename="m_valious";
 m_mat=readmatrix('M_tab.txt');
 save("m_valious","m_mat");

