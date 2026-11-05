a = randi(1e6,1e7,10);
arr = ArrayOfArrays(a);

%%

m = arr.getArray(1:1e6);
m1 = a(1:1e6,:);
e = norm(m1 - m,'f');
M = (reshape(m,10,[]))';
