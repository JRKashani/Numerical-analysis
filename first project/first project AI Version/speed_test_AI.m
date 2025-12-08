x = randi([1 600], 1, 1000)/10;
y = randi(10000, 1, 1000);

tic;
for i = 1:1000
    SolveBesselAI(x(i), y(i));
end
toc;