clear
close all
delete(gcp('nocreate'))

tic
n = 200;
A = 500;
a = zeros(n);
for i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
t = toc;
fprintf('Time with for loop %2.4f \n',t)

tic
parpool("Threads",4)
n = 200;
A = 500;
a = zeros(n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
tpar = toc;
fprintf('Time with parallel for loop %2.4f \n',tpar)