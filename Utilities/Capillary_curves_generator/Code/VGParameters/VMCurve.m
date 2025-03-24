function [interval,Swb,KrwA,KrwB] = VMCurve(kpa,n,beta,lims,pts)
% Function to make the capillary curve and retantion curve
% Ex.:
% [interval,Swb,KrwA,KrwB] = VMCurve(0.5,3,1,[-2,1],1000);


% Find the parameteres of the curves.
m=1-1/n;
pc=1/beta;

% Define the curves.
SeFun = @(p,b,n,m) (1 + (b*p).^n).^(-m);
DSeFun = @(p,b,n,m) ((-m.*n.*p.^(-1)).*(b.*p).^n .* (1+(b.*p).^n).^(-1-m));
DDSeFun = @(p,b,n,m) ((m.*n.*p.^(-2)).*(b.*p).^n .* (1+(b.*p).^n).^(-2-m).*(n.*(m.*(b.*p).^n -1) + (b.*p).^n + 1));

krwFunA = @(p,pc,n,m) (1 + (p/pc).^n).^(-2.5.*m) .* ((1 + (p/pc).^n).^m - ((p/pc).^n).^m).^2;

krwFunB = @(p,beta,kpa,n,m) (1 + (beta*p).^n).^(-kpa.*m) .* (1 - (1 - ((1 + (beta*p).^n).^(-m)).^(1/m)).^m).^2;
% krnFunB = @(p,beta,kpa,n,m) (1-(1 + (beta*p).^n).^(-m)).^kpa .* (1 - (1 + (beta*p).^n).^(-1)).^(2*m);

krwFunC = @(p,pc,kpa,n,m) SeFun(p,pc,n,m).^(kpa) .* (1 - SeFun(p,pc,n,m) .*((1 + (p/pc).^n) -1).^m).^2;
% krwFunD = @(p,beta,kpa,n,m) (1 + (beta*p).^n).^(-kpa.*m) .* (1-((1+(beta*p).^n).^(-m)) .* ((beta*p).^(n*m))).^2;

krwFunD = @(p,b,k,n,m) (1 + (b*p).^n).^(-k*m) .* (1- ((b*p).^(n*m)).*((1+(b*p).^n).^(-m))).^2;
DkrwFunD = @(p,b,k,n,m)  - (m*n./p) .* (((b*p).^(m*n))-(((b*p).^n +1).^m)) .* (((b*p).^n +1).^(-m*(k+2)-1)) ... 
    .*( (k*(b*p).^n).*((b*p).^(m*n) - ((b*p).^n + 1).^m) - 2*((b*p).^(m*n)) );

krwBurdD = @(p,beta,n,m) (1 + (beta*p).^n).^(-2.*m) .* (1- ((beta*p).^(n*m)).*((1+(beta*p).^n).^(-m)));
DkrwBurdD = @(p,beta,n,m) -(m*n*(p.^(-1))) .* ((1+(beta.*p).^n).^(-3*m-1)) ...
    .* (2.*((beta.*p).^n).*((1+(beta.*p).^n).^m) + (beta*p).^(m*n) ...
    - 2.*(beta.*p).^(n*m+n));

% Find the points of the curves.
interval = linspace(10^lims(1),10^lims(2),pts);
Swb = SeFun(interval,beta,n,m);
DSwb = DSeFun(interval,beta,n,m);
DDSwb = DDSeFun(interval,beta,n,m);

dSwb2 = diff(Swb)/((10^lims(2)-10^lims(1))/pts);
ddSwb2 = diff(DSwb)/((10^lims(2)-10^lims(1))/pts);



KrwA = krwFunA(interval,pc,n,m);
KrwB = krwFunB(interval,beta,kpa,n,m);
KrwC = krwFunC(interval,pc,kpa,n,m);

KrwD = krwFunD(interval,beta,kpa,n,m);
dKrwD = DkrwFunD(interval,beta,kpa,n,m);
dKrwD2 = diff(KrwD)/((10^lims(2)-10^lims(1))/pts);

KrwBurdD = krwBurdD(interval,beta,n,m);
dKrwBurdD = DkrwBurdD(interval,beta,n,m);

% KrnB = krnFunB(interval,beta,kpa,n,m);

error_cap = abs(KrwD-KrwBurdD)./KrwD;
error_dk_diff = abs(dKrwD(2:length(interval))-dKrwD2)./abs(dKrwD(2:length(interval)));
error_dS_diff = abs(DSwb(2:length(interval))-dSwb2)./abs(DSwb(2:length(interval)));
error_ddS_diff = abs(DDSwb(2:length(interval))-ddSwb2)./abs(DDSwb(2:length(interval)));

% Plots.
figure();
hold on;
plot(interval,Swb,'r');
title('Capillary curve');
xlabel('Pressure');
ylabel('S_e');
xlim([10^lims(1),10^lims(2)]);
ylim([0, 1]);
hold off;

figure();
hold on;
plot(interval,KrwA,'r');
plot(interval,KrwB,'b');
plot(interval,KrwC,'g');
plot(interval,KrwD,'black');
plot(interval,KrwBurdD,'b--');
title('Relative permeability curve');
xlabel('Pressure');
ylabel('k_r');
xlim([10^lims(1),10^lims(2)]);
ylim([0, 1]);
legend('Reference','MRST','ḾRST2','ḾRST3','Burdine')
hold off;






figure();
hold on;
plot(interval,DSwb,'r');
plot(interval(2:length(interval)),dSwb2,'b');
title('Derivative capillary curve');
xlabel('Pressure');
ylabel('dS_e');
xlim([10^lims(1),10^lims(2)]);
% ylim([0, 1]);
legend('Wolfram','diff')
hold off;


figure();
hold on;
plot(interval,DDSwb,'r');
plot(interval(2:length(interval)),ddSwb2,'b');
title('Second Derivative capillary curve');
xlabel('Pressure');
ylabel('dS_e');
xlim([10^lims(1),10^lims(2)]);
% ylim([0, 1]);
legend('Wolfram','diff')
hold off;


% Plots.
figure();
hold on;
plot(interval,dKrwD,'r');
plot(interval,dKrwBurdD,'b--');
plot(interval(2:length(interval)),dKrwD2,'b');
title('Derivative Relative permeability curve');
xlabel('Pressure');
ylabel('dk_r');
xlim([10^lims(1),10^lims(2)]);
% ylim([0, 1]);
legend('Wolfram','diff','Burdine')
hold off;


figure();
hold on;
plot(interval,error_cap,'red');
plot(interval(2:length(interval)),error_dk_diff,'blue');
plot(interval(2:length(interval)),error_dS_diff,'yellow');
plot(interval(2:length(interval)),error_ddS_diff,'black');
title('Relative Error between the curves');
xlabel('Pressure');
ylabel('error');
xlim([10^lims(1),10^lims(2)]);
legend('k_{error}','dk_{error}','dS_{error}','ddS_{error}')
hold off;

end