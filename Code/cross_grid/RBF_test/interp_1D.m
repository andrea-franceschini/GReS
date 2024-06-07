%% RBF TEST (1d case) discontinous linear shape function integrated within the domain
clear;
close all;

% kink position in the interval -1;1


% gauss points class
% g = Gauss(12,15,1);
% sample points for RBF interpolation
%pts = unique([g.coord' -1 k 1]);
%pts = linspace(-1,1,11);

nInt = 4; % number of interpolation points over each element (excluding nodes)
tol = 1e-5;

%% CASE 1 - interpolate basis function on a single element
a = -1;
b = 1;
%d = (b-a)/1000;
ptsX = linspace(a,b,nInt); 
ptsY = ptsX;
%ptsY = 2*ones(legnth(ptsX))
%ptsY = zeros(length(ptsX),1);
%N = @(x) 1-x.^2;
N = @(x) -0.5*x+0.5;
%N = @(x) sin(x);
vals= N(ptsX);

type = 'gauss';
%% INTERPOLATE USING SELECTED RBF
%r = b-a;
%r = sqrt((max(ptsX)-min(ptsX)).^2 + (max(ptsY)-min(ptsY)).^2);
r = 50;
fiMM = zeros(length(ptsX),length(ptsX));
PM = zeros(length(ptsX),3);
for i = 1:length(ptsX)
    dist = (ptsX - ptsX(i)).^2 + (ptsY - ptsY(i)).^2;
    dist = sqrt(dist);
    rbf = computeRBFentries(dist,type,r);
    fiMM(:,i) = rbf;
end

% polynomial term
for i = 1:length(ptsX)
    PM(i,:) = [1 ptsX(i) ptsY(i)];
end

mat = [fiMM PM;
    PM' zeros(3,3)];
w = mat\[vals';zeros(3,1)];
wf = w(1:length(ptsX));
alpha = w(length(ptsX)+1:end);
%w1 = fiMM\ones(length(ptsX),1);
% sample2 = linspace(-1,1,100);

% compute normal 
n = [1,-1];
k = 0.5;

% curved slave sample to check the effect of non aligned mesh
sampX = linspace(-1.2,1.2,20);
sampleX = sampX+k*n(1);
sampleY = sampX + k*n(2);
iiVec = []; jjVec =[]; rbfVec = [];
fiNM = zeros(length(sampX),length(ptsX));
PN = zeros(length(sampX),3);
for i = 1:length(sampleX)
    dist = sqrt((ptsX - sampleX(i)).^2 + (ptsY - sampleY(i)).^2);
    rbf = computeRBFentries(dist,type,r);
    fiNM(i,:) = rbf;
end
for i = 1:length(ptsX)
    PM(i,:) = [1 sampleX(i) sampleY(i)];
end

mat2 = [fiNM PN];

%fiNM = sparse(iiVec,jjVec,rbfVec,length(sampleX),length(ptsX));
%valsOut = (fiNM*wf)./(fiNM*w1);
vals = (mat2*w);
valsOut = vals(1:end-2);
%vals2_noRL = (fiNM*wf);
%pts_c1 = pts;




%%

% %% CASE 2 - interpolate basis function over the support of a node
% 
% a = -1;
% b = 1;
% k = 0;
% pts1 = [a linspace(a+0.01,k-0.01,nInt) k];
% pts2 = [k linspace(k+0.01,b-0.01,nInt) b];
% pts = unique([pts1 pts2]);
% % pts = [linspace(-2,a-0.01,nInt)]
% N1 = @(x) 1/(k-a)*x - a/(k-a);
% N2 = @(x) 1/(k-b)*x - b/(k-b);
% vals1 = N1(pts1);
% vals2 = N2(pts2);
% vals = [vals1(1:end) vals2(2:end)];
% 
% r = b-a;
% iiVec = []; jjVec =[]; rbfVec = [];
% for i = 1:length(pts)
%     dist = pts - pts(i);
%     dist = sqrt(dist.^2);
%     rbf = pos(1-dist./r).^4.*(1+4*dist./r);
%     iiVec = [iiVec; repmat(i,length(pts),1)];
%     jjVec = [jjVec; (1:length(pts))'];
%     rbfVec = [rbfVec; rbf'];
% end
% fiMM = sparse(iiVec,jjVec,rbfVec,length(pts),length(pts));
% 
% wf = fiMM\vals';
% w1 = fiMM\ones(length(pts),1);
% sample2 = linspace(-1.3,1.3,200); 
% iiVec = []; jjVec =[]; rbfVec = [];
% for i = 1:length(sample2)
%     dist = pts - sample2(i);
%     dist = sqrt(dist.^2);
%     rbf = pos(1-dist./r).^4.*(1+4*dist./r);
%     iiVec = [iiVec; repmat(i,length(pts),1)];
%     jjVec = [jjVec; (1:length(pts))'];
%     rbfVec = [rbfVec; rbf'];
% end
% fiNM = sparse(iiVec,jjVec,rbfVec,length(sample2),length(pts));
% vals2_c2 = (fiNM*wf)./(fiNM*w1);
% pts_c2 = pts;
% sample_c2 = sample2;
% 
% 
% %% CASE 3 - interpolate basis function over the support of a node and all its neighbors
% a = -1;
% b = 1;
% k = 0;
% pts1 = [a linspace(a+0.01,k-0.01,nInt) k];WENDLAND C2 FUNCTIONS
% pts2 = [k linspace(k+0.01,b-0.01,nInt) b];
% pts = unique([pts1 pts2]);
% pts = [linspace(a-1,a-0.01,nInt) pts linspace(b+0.01,b+1,nInt)];
% % pts = [linspace(-2,a-0.01,nInt)]
% N1 = @(x) 1/(k-a)*x - a/(k-a);
% N2 = @(x) 1/(k-b)*x - b/(k-b);
% vals1 = N1(pts1);
% vals2 = N2(pts2);
% vals = [zeros(1,nInt) vals1(1:end) vals2(2:end) zeros(1,nInt)];
% r = b+2;
% iiVec = []; jjVec =[]; rbfVec = [];
% for i = 1:length(pts)
%     dist = pts - pts(i);
%     dist = sqrt(dist.^2);
%     rbf = pos(1-dist./r).^4.*(1+4*dist./r);
%     iiVec = [iiVec; repmat(i,length(pts),1)];
%     jjVec = [jjVec; (1:length(pts))'];
%     rbfVec = [rbfVec; rbf'];
% end
% fiMM = sparse(iiVec,jjVec,rbfVec,length(pts),length(pts));
% 
% wf = fiMM\vals';
% w1 = fiMM\ones(length(pts),1);
% sample2 = linspace(-2,2,200); 
% iiVec = []; jjVec =[]; rbfVec = [];
% for i = 1:length(sample2)
%     dist = pts - sample2(i);
%     dist = sqrt(dist.^2);
%     rbf = pos(1-dist./r).^4.*(1+4*dist./r);
%     iiVec = [iiVec; repmat(i,length(pts),1)];
%     jjVec = [jjVec; (1:length(pts))'];
%     rbfVec = [rbfVec; rbf'];
% end
% fiNM = sparse(iiVec,jjVec,rbfVec,length(sample2),length(pts));
% vals2_c3 = (fiNM*wf)./(fiNM*w1);
% pts_c3 = pts;
% sample_c3 = sample2;


%% plot 
xval = linspace(-1,1,100);
yval = N(xval);
%plot(ptsX, vals, 'o', 'LineWidth',1)
plot(xval,yval,'k-')
hold on
plot(sampleX,valsOut, 'r*')
%plot(sample_c1,vals2_noRL, 'g*')
% plot(sample_c2,vals2_c2, 'b^')
% plot(sample_c3,vals2_c3, 'gs')
y = N(sampX');

%% measure L2 error
L2 = norm(2/length(xval)*(valsOut - y),2);

%% error comparison 

% h = 2/(length(sample2)-1);
% 
% valsEx1 = N1(sample_c1);
% err1 = (valsEx1'-vals_c1).^2;
% err1 = sqrt(sum(err1*h));
% 
% vEx1 = N1(sample_c2);
% vEx2 = N2(sample_c2);
% err2 = (valsEx'-valsNoScale).^2;
% err2 = sqrt(sum(errNoScale*h));


%% INTEGRATION TEST
% N = @(x) 0.5*x+0.5;
% int = [];
% ex = (int(2)-int(1))*N(mean(int));
% 
% 
% weights = sort(weights);
% if rem(length(weights),2)~=0
%     weights = [weights flip(weights(1:end-1))];
% else
%     weights = [weights flip(weights(1:end))];    
% end
% gauss_int = 0;
% for i = 1:length(gPoints)
%     if int(1)<gPoints(i) && gPoints(i)<int(2)
%         gauss_int = gauss_int + N(gPoints(i))*weights(i);
%     end
% end
% err = ex-gauss_int;

