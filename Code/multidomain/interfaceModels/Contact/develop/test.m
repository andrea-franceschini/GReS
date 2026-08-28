v1 = randi(10,30,1);
v2 = randi(7,20,1);

for i = 1:1e6
  v1u = unique(v1);
  v2u = unique(v2);
  v3 = setdiff(v1u,v2u);
end