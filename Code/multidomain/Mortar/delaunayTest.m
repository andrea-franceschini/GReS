rng default;
x = rand([20,1]);
y = rand([20,1]);
DT = delaunay(x,y);
triplot(DT,x,y);