%% RBF TEST (1d case) discontinous linear shape function integrated within the domain
clear;
close all;

% kink position in the interval -1;1
k = 0;

% gauss points class
g = Gauss(12,15,1);
% sample points for RBF interpolation
%pts = unique([g.coord' -1 k 1]);
%pts = linspace(-1,1,11);

nInts = 5;
pts = [-1 -0.05 0 0.05 1];
pts = [linspace(pts(1), pts(2),nInts), pts(3), linspace(pts(4),pts(5),nInts)];
% piecewise linear functions
c1 = 1/(1+k);
c2 = 1/(1-k);
N1 = @(x) c1*x + c1;
N2 = @(x) -c2*x + c2;
vals = zeros(length(pts),1);
for i = 1:length(pts)
    if -1.01 < pts(i) && pts(i) < k
        % vals(i) = N1(sample(i));
        vals(i) = N1(pts(i));
    elseif k-0.001 < pts(i) && pts(i) < k+0.001
        % vals(i) = 1;
        vals(i) = 1;
    elseif k < pts(i) && pts(i) < 1.01
         vals(i) = N2(pts(i));
    end
end
r = 2;
iiVec = []; jjVec =[]; rbfVec = [];
tic
for i = 1:length(pts)
    dist = pts - pts(i);
    dist = sqrt(dist.^2);
    rbf = pos(1-dist./r).^4.*(1+4*dist./r);
    iiVec = [iiVec; repmat(i,length(pts),1)];
    jjVec = [jjVec; (1:length(pts))'];
    rbfVec = [rbfVec; rbf'];
end
fiMM = sparse(iiVec,jjVec,rbfVec,length(pts),length(pts));


wf = fiMM\vals;
w1 = fiMM\ones(length(pts),1);
toc
sample2 = linspace(-5,5,500); 
iiVec = []; jjVec =[]; rbfVec = [];
for i = 1:length(sample2)
    dist = pts - sample2(i);
    dist = sqrt(dist.^2);
    rbf = pos(1-dist./r).^4.*(1+4*dist./r);
    iiVec = [iiVec; repmat(i,length(pts),1)];
    jjVec = [jjVec; (1:length(pts))'];
    rbfVec = [rbfVec; rbf'];
end
fiNM = sparse(iiVec,jjVec,rbfVec,length(sample2),length(pts));
vals2 = (fiNM*wf)./(fiNM*w1);
plot(pts, vals, 'k-o')
hold on
plot(sample2,vals2, 'r')
xlim([-5 5])
ylim([-0.6 1.1])


%% INTEGRATION TEST
N = @(x) 0.5*x+0.5;
int = [];
ex = (int(2)-int(1))*N(mean(int));


weights = sort(weights);
if rem(length(weights),2)~=0
    weights = [weights flip(weights(1:end-1))];
else
    weights = [weights flip(weights(1:end))];    
end
gauss_int = 0;
for i = 1:length(gPoints)
    if int(1)<gPoints(i) && gPoints(i)<int(2)
        gauss_int = gauss_int + N(gPoints(i))*weights(i);
    end
end
err = ex-gauss_int;
