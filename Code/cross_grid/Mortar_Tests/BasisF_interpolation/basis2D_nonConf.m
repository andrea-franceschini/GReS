% study the behaviour of different RBFs when interpolating a basis function
% defined over the interval -1,1.
% we use uniformly distributed interpolation points (also Gauss point
% position can be an option)
% we check the RMS error for fixed r and different N points and fixed N and
% different r
% we study the behaviour for different families of RBF
% we come up with the plot of chapter 3 of the paper
close all; clear
warning('off','MATLAB:nearlySingularMatrix');
a = -1; b = 1;
N = 2;
%f = @(x) 0.5*(1-x(:,1).^2).*(1+x(:,2));
%f = @(x) 0.25*(1+x(:,1)).*(1+x(:,2)) + 0*x(:,3);
%f = @(x) 1-max(abs(x(:,1)),abs(x(:,2)));
f = @(x) 0.5-0.5*x(:,1);
type = 'gauss';
fac = 1;
intPts = linspace(a,b,N);
%intPts = sin(pi*intPts/2);
[y, x] = meshgrid(intPts, intPts);
z = zeros(length(x));
intPts = [x(:), y(:), z(:)];
vals = f(intPts);
intPts(:,2) = fac*intPts(:,2);
r = sqrt(1+fac^2)*(b-a);
fiMM = zeros(length(intPts),length(intPts));
for ii = 1:length(intPts)
    dist = (intPts(:,1) - intPts(ii,1)).^2 + (intPts(:,2) - intPts(ii,2)).^2 + (intPts(:,3) - intPts(ii,3)).^2;
    dist = sqrt(dist);
    fiMM(ii,:) = computeRBFentries(dist,type,r);
end
wf = fiMM\vals;
w1 = fiMM\ones(length(intPts),1);

N  = 100;
samp = linspace(-50,50,N);
[ys, xs] = meshgrid(samp,samp);
zs = 0*xs;
samp = [xs(:), ys(:), zs(:)];
%samp(:,2) = fac*samp(:,2);
fiNM = zeros(length(samp),length(intPts));
for ii = 1:length(samp)
    dist = (intPts(:,1) - samp(ii,1)).^2 + (intPts(:,2) - samp(ii,2)).^2 + (intPts(:,3) - samp(ii,3)).^2;
    dist = sqrt(dist);
    fiNM(ii,:) = computeRBFentries(dist,type,r);
end
valsOut = (fiNM*wf)./(fiNM*w1);
%valsOut = (fiNM*wf);
yval = f(samp);
% error
L2 = norm((1/length(samp))*(valsOut - yval),2);

scatter3(samp(:,1),samp(:,2),valsOut)
hold on
surf(xs,ys,zeros(length(xs),length(ys)))


figure(2)
surf(x,y,z)
hold on
surf(xs,ys,zs)




