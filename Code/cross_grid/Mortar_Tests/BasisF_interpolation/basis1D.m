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
a = -1;
b = 1;
f = @(x) 0.5*x.*(1+x);
%f = @(x) -0.5 + 0.5*x;
N = 4:2:24;
%fac = [2];
%type = 'gauss';
typeStr = ["gauss","imq","wendland"];
L2 = zeros(numel(N),numel(typeStr));
n_cond = zeros(numel(N),numel(typeStr));
mark = {'o','*','^'};
k = 0;
%for j = 1:length(fac)
for type = typeStr
    k = k+1;
    for i = 1:length(N)
        xInt = linspace(a,b,N(i));
        vals= f(xInt);
        r = (b-a);
        fiMM = zeros(length(xInt),length(xInt));
        for ii = 1:length(xInt)
            dist = (xInt - xInt(ii)).^2;
            dist = sqrt(dist);
            fiMM(ii,:) = computeRBFentries(dist,type,r);
        end
        wf = fiMM\vals';
        w1 = fiMM\ones(length(xInt),1);

        sampX = linspace(a,b,40);
        fiNM = zeros(length(sampX),length(xInt));
        for ii = 1:length(sampX)
            dist = (xInt - sampX(ii)).^2;
            dist = sqrt(dist);
            fiNM(ii,:) = computeRBFentries(dist,type,r);
        end
        valsOut = (fiNM*wf)./(fiNM*w1);
        %valsOut = fiNM*wf;
        y = f(sampX');
        % error
        L2(i,k) = norm((1/length(sampX))*(valsOut - y),2);
        n_cond(i,k) = cond(fiMM);
    end
    % figure(1)
    % semilogy(N,L2,'-k','LineWidth',1,'Marker',mark{k})
    % xlabel('M')
    % ylabel('RMSE')
    % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    % ax = get(gca,'XTickLabel');
    % set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
    % grid on
    % hold on
    % figure(2)
    % semilogy(N,n_cond,'-k','LineWidth',1,'Marker',mark{k})
    % xlabel('M')
    % ylabel('Condition number')
    % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    % ax = get(gca,'XTickLabel');
    % set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
    % grid on
    % hold on

    fid = fopen("c1_L2",'w');
    fprintf(fid,'%2.6e %2.6e %2.6e\n',L2');
    fid = fopen("c1_cond.dat",'w');
    fprintf(fid,'%2.6e %2.6e %2.6e\n',n_cond');
end
%legend(typeStr)

%%
% x = linspace(-1,1,50);
% bf = @(x) 0.5*x.*(1+x);
% plot(x,bf(x),'b-','LineWidth',1.5);
% %colormap
% xs = linspace(-1,1,6);
% hold on
% scatter(xs,zeros(length(xs),1),"red")
% scatter(xs,bf(xs),"filled","red")
% xlim([-1.2 1.2])
% ylim([-0.2 1.2])
% set(gca,'XTick',[])
% set(gca,'YTick',[])













