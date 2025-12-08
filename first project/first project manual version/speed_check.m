x = randi([1 600], 1, 1000)/10;
y = randi(10000, 1, 1000);

tic;
for i = 1:1000
    project1(x(i), y(i));
end
toc;