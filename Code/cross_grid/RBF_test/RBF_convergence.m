%% RBF TEST: accuracy and condition number vs. radius
clear;
close all;

% kink position in the interval -1;1

%% CASE 1 - interpolate basis function on a single element
a = -1;
b = 1;
n_refs = 7;
L2 = zeros(n_refs,1);
n_cond = zeros(n_refs,1);
nInts = [4 8 16 32 64 128 256];
%d = (b-a)/1000;
for ii=1:length(nInts)
    nInt = nInts(ii);
    ptsX = linspace(a,b,nInt);
    ptsY = 0*ptsX;
    %ptsY = 2*ones(legnth(ptsX))
    %ptsY = zeros(length(ptsX),1);
    %N = @(x) -0.5*x+0.5;
    N = @(x) 0.5*x.*(1+x);
    %N = @(x) sin(x);
    vals= N(ptsX);
    %vals = ones(1,length(ptsX));
    type = 'gauss';
    %% INTERPOLATE USING SELECTED RBF
    %r = b-a;
    %r = sqrt((max(ptsX)-min(ptsX)).^2 + (max(ptsY)-min(ptsY)).^2);
    %rad = (2/length(ptsX))*linspace(0.5,20,100);
    %n_cond = zeros(length(rad),1);
    r = b-a;
    iiVec = []; jjVec =[]; rbfVec = [];
    for i = 1:length(ptsX)
        dist = (ptsX - ptsX(i)).^2 + (ptsY - ptsY(i)).^2;
        dist = sqrt(dist);
        rbf = computeRBFentries(dist,type,r);
        iiVec = [iiVec; repmat(i,length(ptsX),1)];
        jjVec = [jjVec; (1:length(ptsX))'];
        rbfVec = [rbfVec; rbf'];
    end
    fiMM = sparse(iiVec,jjVec,rbfVec,length(ptsX),length(ptsX));

    wf = fiMM\vals';
    w1 = fiMM\ones(length(ptsX),1);
    % sample2 = linspace(-1,1,100);

    % compute normal
    n = [1,-1];
    k = 0.5;

    % curved slave sample to check the effect of non aligned mesh
    sampX = linspace(a,b,200);
    sampleX = sampX;
    sampleY = 0*sampX;
    iiVec = []; jjVec =[]; rbfVec = [];
    for i = 1:length(sampleX)
        dist = sqrt((ptsX - sampleX(i)).^2 + (ptsY - sampleY(i)).^2);
        rbf = computeRBFentries(dist,type,r);
        iiVec = [iiVec; repmat(i,length(ptsX),1)];
        jjVec = [jjVec; (1:length(ptsX))'];
        rbfVec = [rbfVec; rbf'];
    end
    fiNM = sparse(iiVec,jjVec,rbfVec,length(sampleX),length(ptsX));
    valsOut = (fiNM*wf)./(fiNM*w1);
    %valsOut = (fiNM*wf);
    %vals2_noRL = (fiNM*wf);
    %pts_c1 = pts;
    y = N(sampX');

    %% measure L2 error
    L2(ii) = norm(2/length(sampX)*(valsOut - y),2);
    n_cond(ii) = cond(fiMM);
end
figure(1)
semilogy(4*2.^(1:n_refs),L2,'k-o', 'LineWidth',1)
xlabel('nPoints')
ylabel('L2 approximation error')

figure(2)
semilogy(4*2.^(1:n_refs),n_cond,'k-o','LineWidth',1)
xlabel('nPoints')
ylabel('Condtion numb. interpolation mat')

%%
g = load("gauss.mat");
w = load("wendland.mat");
L2_g = g.L2;
L2_w = w.L2;
cond_g = g.n_cond;
cond_w = w.n_cond;

figure(1)
loglog(4*2.^(1:n_refs),L2_g,'r-o','LineWidth',1)
hold on
loglog(4*2.^(1:n_refs),L2_w,'b-o','LineWidth',1)
xlabel('nPoints')
ylabel('L2 interp. error')
legend('gauss','wendland')

figure(2)
loglog(4*2.^(1:n_refs),cond_g,'r-o','LineWidth',1)
hold on
loglog(4*2.^(1:n_refs),cond_w,'b-o','LineWidth',1)
xlabel('nPoints')
ylabel('condition number')
legend('gauss','wendland')









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
% xval = linspace(-1,1,100);
% yval = N(xval);
% %plot(ptsX, vals, 'o', 'LineWidth',1)
% plot(xval,yval,'k-')
% hold on
% plot(sampleX,valsOut, 'r*')
%plot(sample_c1,vals2_noRL, 'g*')
% plot(sample_c2,vals2_c2, 'b^')
% plot(sample_c3,vals2_c3, 'gs')


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

