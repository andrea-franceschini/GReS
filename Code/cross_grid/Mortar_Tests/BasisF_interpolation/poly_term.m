% study the behaviour of different RBFs when interpolating a basis function
% defined over the interval -1,1.
% we use uniformly distributed interpolation points (also Gauss point
% position can be an option)
% we check the RMS error for fixed r and different N points and fixed N and
% different r
% we study the behaviour for different families of RBF
% we come up with the plot of chapter 3 of the paper

% INFLUENCE OF THE POLYNOMIAL TERM
close all; clear
warning('off','MATLAB:nearlySingularMatrix');
a = -3;
b = 3;
%f = @(x) 1+0*x;
%f = @(x) 0.5*x.*(1+x);
f = @(x) -0.5 + 0.5*x;
N = 5;
%fac = [2];
%L2 = zeros(numel(N),1);
% n_cond = zeros(numel(N),1);
type = 'gauss';

pts = linspace(-1,1,N);
int = ref2nod(pts,[-3,-3],[3,3]);
int = int';
xInt = int(1,:);
yInt = int(2,:);
vals= f(pts');
r = (b-a);
fiMM = zeros(length(xInt),length(xInt));
for ii = 1:length(xInt)
    dist = (xInt - xInt(ii)).^2 + (yInt - yInt(ii)).^2;
    dist = sqrt(dist);
    fiMM(ii,:) = computeRBFentries(dist,type,r);
end
P = [ones(length(xInt),1) pts'];
wf = [fiMM P;P' zeros(2)]\[vals; zeros(2,1)];
%w1 = [fiMM P;P' zeros(3)]\[ones(length(xInt),1); zeros(3,1)];

samp = linspace(-1,1,50);
int = ref2nod(samp,[-3.2,-3.2],[3.2,3.2]);
int = int';
sampX = int(1,:)+1;
sampY = int(2,:)-1;
fiNM = zeros(length(samp),length(pts));
for ii = 1:length(sampX)
    dist = (xInt - sampX(ii)).^2 + (yInt - sampY(ii)).^2;
    dist = sqrt(dist);
    fiNM(ii,:) = computeRBFentries(dist,type,r);
end
Pn = [ones(length(samp),1) samp'];
A = [fiNM Pn];
%valsOut = (A*wf)./(A*w1);
valsOut = A*wf;
y = f(sampX');
% error
L2 = norm((1/length(sampX))*(valsOut - y),2);
n_cond = cond(fiMM);
%xP = linspace(-1,1,500);

plot(sampX,valsOut, 'b*')
hold on
%plot(xP,f(xP), 'r')









