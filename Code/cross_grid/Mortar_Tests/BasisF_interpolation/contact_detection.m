pts = linspace(-2,2,50);
linspace(-2,2,50);
f = @(x,y) 1-x;
[x,y] = meshgrid(pts,pts);

surf(x,y,f(x,y))