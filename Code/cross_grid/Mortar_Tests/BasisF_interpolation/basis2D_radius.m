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
N = 4:2:24;
L2 = zeros(numel(N),1);
%n_cond = zeros(numel(N),1);
f = @(x) 0.5*(1-x(:,1).^2).*(1+x(:,2));
%f = @(x) 0.25*(1+x(:,1)).*(1+x(:,2));
%for j = 1:length(fac)
fac = 1;
mark = {'o','*','^'};
k = 10;
type = 'gauss';
for i = 1:length(N)
    intPts = linspace(a,b, N(i));
    [y, x] = meshgrid(intPts, intPts);
    intPts = [x(:), y(:)];
    vals= f(intPts);
    % FIXED RADIUS
    rfix = sqrt(1+fac^2)*(b-a);
    % distance matrix
    D = zeros(length(intPts),length(intPts)-1);
    % ADAPTIVE RADIUS BASED ON THE FILL DISTANCE
    for ii = 1:length(intPts)
        dist = (intPts([1:ii-1 ii+1:end],1) - intPts(ii,1)).^2 + (intPts([1:ii-1 ii+1:end],2) - intPts(ii,2)).^2;
        dist = sqrt(dist);
        D(ii,:) = dist;
    end
    r = max(min(D));
    r = k*r;
    fiMM = zeros(length(intPts),length(intPts));
    for ii = 1:length(intPts)
        dist = (intPts(:,1) - intPts(ii,1)).^2 + (intPts(:,2) - intPts(ii,2)).^2;
        dist = sqrt(dist);
        fiMM(ii,:) = computeRBFentries(dist,type,r);
    end
    wf = fiMM\vals;
    w1 = fiMM\ones(length(intPts),1);
    h = haltonset(2,'Skip',1e3,'Leap',1e2); % ad-hoc interpolation
    samp = net(h,50);
    [y, x] = meshgrid(samp(:,1),samp(:,2));
    samp = [x(:), y(:)];
    sampUndef = samp;
    fiNM = zeros(length(samp),length(intPts));
    for ii = 1:length(samp)
        dist = (intPts(:,1) - samp(ii,1)).^2 + (intPts(:,2) - samp(ii,2)).^2;
        dist = sqrt(dist);
        fiNM(ii,:) = computeRBFentries(dist,type,r);
    end
    valsOut = (fiNM*wf)./(fiNM*w1);
    %valsOut = (fiNM*wf);
    y = f(sampUndef);
    % error
    L2(i) = norm((1/length(samp))*(valsOut - y),2);
    %n_cond(i) = cond(fiMM);
end
% set(axL2, 'YScale', 'log')
% set(axL2, 'XScale', 'log')
% set(axH1, 'YScale', 'log')
% set(axH1, 'XScale', 'log')

switch type
    case 'wendland'
        name = strcat('graph_rad/fill',num2str(k),'W.dat');
        fid = fopen(name,'w');
        fprintf(fid,'%2.6e \n',L2);
    case 'gauss'
        name = strcat('graph_rad/fill',num2str(k),'G.dat');
        fid = fopen(name,'w');
        fprintf(fid,'%2.6e \n',L2);
end

%%










% 
%     figure(1)
%     semilogy(N,L2,'k-','LineWidth',1,'Marker',mark{c})
%     xlabel('M')
%     ylabel('RMSE')
%     legend('Uniform grid', 'Modified')
%     %legend('Gaussian Splines', 'IMQ', 'Wendland')
%     set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
%     ax = get(gca,'XTickLabel');
%     %set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%     axis tight
%     grid on
%     hold on
%     figure(2)
%     semilogy(N,n_cond,'k-','LineWidth',1,'Marker',mark{c})
%     xlabel('M')
%     ylabel('Condition number')
%     %legend('Gaussian Splines', 'IMQ', 'Wendland')
%     set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
%     ax = get(gca,'XTickLabel');
%     %set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 10)
%     axis tight
%     grid on
%     hold on
% end
% %end

    % 
    % for i = 1:length(N)
    %     intPts = linspace(a,b, N(i));
    %     %intPts = sin(pi*intPts/2);
    %     [y, x] = meshgrid(intPts, intPts);
    %     intPts = [x(:), y(:)];
    %     vals= f(intPts);
    %     r = sqrt(2)*(b-a);
    %     fiMM = zeros(length(intPts),length(intPts));
    %     for ii = 1:length(intPts)
    %         dist = (intPts(:,1) - intPts(ii,1)).^2 + (intPts(:,2) - intPts(ii,2)).^2;
    %         dist = sqrt(dist);
    %         fiMM(ii,:) = computeRBFentries(dist,type,r);
    %     end
    %     wf = fiMM\vals;
    %     w1 = fiMM\ones(length(intPts),1);
    % 
    %     samp = [-1 1]; % ad-hoc interpolation
    %     samp = linspace(samp(1), samp(2), 40);
    %     [y, x] = meshgrid(samp,samp);
    %     samp = [x(:), y(:)];
    %     fiNM = zeros(length(samp),length(intPts));
    %     for ii = 1:length(samp)
    %         dist = (intPts(:,1) - samp(ii,1)).^2 + (intPts(:,2) - samp(ii,2)).^2;
    %         dist = sqrt(dist);
    %         fiNM(ii,:) = computeRBFentries(dist,type,r);
    %     end
    %     valsOut = (fiNM*wf)./(fiNM*w1);
    %     %valsOut = (fiNM*wf);
    %     y = f(samp);
    %     % error
    %     L2(i) = norm((1/length(samp))*(valsOut - y),2);
    %     n_cond(i) = cond(fiMM);
    % end
    % figure(1)
    % semilogy(N,L2,'-rs','LineWidth',1)
    % xlabel('M')
    % ylabel('RMSE')
    % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    % ax = get(gca,'XTickLabel');
    % set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 10)
    % hold on
    % figure(2)
    % semilogy(N,n_cond,'-rs','LineWidth',1)
    % xlabel('M')
    % ylabel('Condition number')
    % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    % ax = get(gca,'XTickLabel');
    % set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 10)
    % hold on

    % %% plotting 2D basis function in the master elements
    % figure(1)
    % [X,Y] = meshgrid(-1:.01:1);
    % bf = @(x,y) 0.5*(1-x.^2).*(1+y);
    % Z = bf(X,Y);
    % s = surf(X,Y,Z);
    % s.EdgeColor = 'none';
    % s.FaceLighting = "gouraud";
    % hold on
    % [xs,ys] = meshgrid(linspace(-1,1,6));
    % scatter3(xs,ys,-0.5*ones(size(xs,1)),"red")
    % scatter3(xs,ys,bf(xs,ys),"filled","red")
    % %scatter3(xs(:),ys(:),zeros(length(xs(:),1)),"filled")
    % xlim([-2 1.5])
    % ylim([-2 1.5])
    % set(gca,'XTick',[])
    % set(gca,'YTick',[])
    % set(gca,'ZTick',[])
    % 
    % 








