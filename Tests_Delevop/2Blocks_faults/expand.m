function Aout = expand(A,n)
% expand input matrix by n times in each direction
nr = size(A,1);
nc = size(A,2);

Aout = zeros(n*nr,n*nc);

for i = 1:n
   Aout(i:n:end,i:n:end) = A;
end

end

