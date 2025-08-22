x1 = -2;
x2 = 5;
dual1 = @(x) 0.5*(1-3*x);
dual2 = @(x) 0.5*(1+3*x);

xVal = linspace(-1,1,100);
yVal = dual1(xVal)*x1+dual2(xVal)*x2;

