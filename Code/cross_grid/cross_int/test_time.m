tic
ns = 1000;
for i = 1:ns
    A = rand(32);
    b = ones(32,1);
    c = A*b;
end
t = toc;