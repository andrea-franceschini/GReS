% list of random unique entities
N = 100;
v = (randperm(5 * N, N))';
id = rand(10*N,1) > 0.5;

d1 = zeros(numel(v),1);
d2 = zeros(numel(v),1);


% solution 1: sorting
tic
[vs,ids] = sort(v); % sort array to speed up
s = 0;
k = 0;
iold = 0;
for j = ids'
   s = s+sum(id(iold+1:vs(k+1)));
   iold = vs(k+1);
   d2(j) = s;
   k = k+1;
end
t1 = toc;

% solution 2: no sorting
tic
for k = 1:numel(v)
   d1(k) = sum(id(1:v(k)));
end
t2 = toc;
