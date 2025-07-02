function varargout = quadplot(Q,varargin)

%% QUADPLOT Plots the points of a quadrature rule 
%    QUADPLOT(Q) plots the points of quadrature rule Q, which is a 
%    structure with fields Q.Points, Q.Weights and Q.Properties (see, for 
%    example, quadGaussLegendre and quadtriangle). 
%
%    Q = QUADPLOT(d,Name,Value) same as above but with additional control 
%    over the plot features. Options include specification of the following 
%    Name-Value pair arguments:
%    ______________________________________________________________________
%    'PlotTitle' -- Controls plot title
%    'on' (default) | 'off' 
%    ----------------------------------------------------------------------
%    Turns 'on' (default) or 'off' the plot title. The default title plot
%    includes the number of points, the type, and degree of the quadrature
%    rule.
%    ______________________________________________________________________
%    'PointLabels' -- Controls point labels
%    'off' (default) | 'on' 
%    ----------------------------------------------------------------------  
%    Turns 'off' (default) or 'on' point labels that indicate the index of 
%    the point in Q.Points, i.e., coordinate(s) of the point labeled i are 
%    Q.Points(i,:).
%
%    See also quadGaussJacobi, quadGaussLobatto, quadsquare, quadtriangle
%

%% Validate input

vQ = @(x)validateattributes(x,{'struct'},{});
vP = @(x)validateattributes(x,{'char'},{});
ip = inputParser;
ip.addRequired('Q',vQ);
ip.addParameter('PointLabels','off',vP);
ip.addParameter('PlotTitle','on',vP);
ip.parse(Q,varargin{:});
ip.Results; 
PointLabels = ip.Results.PointLabels; PlotTitle = ip.Results.PlotTitle;

if all(isfield(Q,{'Points','Weights','Properties'}))
    if ~all(isfield(Q.Properties,{'Degree','Type','Domain'}))
        error('Q.Properties must have fields .Degree, .Type, and .Domain');
    end
else
    error(['Input Q must be a structure with fields .Points, .Weights,',...
        ' and .Properties']);
end

%% Plot quadrature rule
n = size(Q.Points,1);
% Set title
if  strcmpi(PlotTitle,'on')
    title({[num2str(n),'-Point ',Q.Properties.Type,' quadrature rule'],...
            ['(Degree ',num2str(Q.Properties.Degree),')'],'',''})
end
hold on; box on; grid on; 
myblue = [0.2422 0.1504 0.6603];
switch size(Q.Points,2)
    case 1
        %% 1D plot    
        a = Q.Properties.Domain(1); b = Q.Properties.Domain(2);
        text(a,-(b-a)/6.5,num2str(a),'HorizontalAlignment','center');
        text(b,-(b-a)/6.5,num2str(b),'HorizontalAlignment','center');
        text((a+b)/2,-(b-a)/6.5,num2str((a+b)/2),'HorizontalAlignment','center');
        plot([a, b],[0 0],'k-','LineWidth',2)
        plot([a, a],[-(b-a)/40 (b-a)/40],'k-','LineWidth',2)
        plot([b, b],[-(b-a)/40 (b-a)/40],'k-','LineWidth',2)
        plot([b, b],[-(b-a)/8 (b-a)/8],'k:')
        plot([a, a],[-(b-a)/8 (b-a)/8],'k:')
        plot([(a+b)/2,(a+b)/2],[-(b-a)/8 (b-a)/8],'k:')
        plot(Q.Points,zeros(size(Q.Points,1)),'Marker','o',...
            'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w',...
            'LineStyle','none');
        plot(Q.Points,zeros(size(Q.Points,1)),'Marker','o',...
            'MarkerSize',4.5,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue,...
            'LineStyle','none'); 
        ax = gca;
        ax.Color = 'none';
        ax.XAxis.Color = 'black';
        ax.XAxisLocation = 'top';
        ax.TickDir = 'both';
        ax.XAxis.TickValues = Q.Points;
        ax.YAxis.Color = 'none';
        if strcmpi(PointLabels,'on')
            ax.XTickLabel = cellstr(strsplit(num2str(1:n)));
        else
            ax.XTickLabel = [];
        end
        daspect([1 1 1]);
        axis([a b -(b-a)/8 (b-a)/8]);
    case 2
        %% 2D plot
        domain.FaceColor = [0.831,0.816,0.784];
        domain.EdgeColor = 'k';
        domain.LineWidth = 3;
        if isempty(Q.Properties.Domain)
            domain.vertices = [ -1 -1; 1 -1; 0 -1+sqrt(3) ];
            Qb = [ Q.Points(:,1), Q.Points(:,2),1-Q.Points(:,1)-Q.Points(:,2)];
            Q.Points = barycentricToCartesian(triangulation([2 3 1],...
               domain.vertices),ones(n,1),Qb);
        else
            domain.vertices = Q.Properties.Domain;
        end
        a = domain.vertices(1,:); b = domain.vertices(2,:); c = domain.vertices(3,:); 
        domain.faces = 1:size(domain.vertices,1);
        patch(domain);
        ax = gca;
        ax.Color = 'none';
        hold on
        switch size(domain.vertices,1)
            case 3 % Triangle
                switch Q.Properties.Type
                    case 'product'
                        j = 1; 
                        m = sqrt(n);
                        for i = 1:m
                            plot([Q.Points(j,1),Q.Points(j+m-1)],...
                                [Q.Points(j,2),Q.Points(j+m-1,2)],'Color',...
                                [0.65 0.65 0.65],'LineWidth',0.1)
                            plot([Q.Points(i,1),Q.Points(n-m+i,1)],...
                                [Q.Points(i,2),Q.Points(n-m+i,2)],'Color',...
                                [0.65 0.65 0.65],'LineWidth',0.1)
                            j = j+m;                           
                        end
                    case 'nonproduct'
                        % Plot lines related to symmetry
                        plot([a(1),(b(1)+c(1))/2],[a(2),(b(2)+c(2))/2],...
                            'Color',[0.65 0.65 0.65],'LineWidth',0.1);
                        plot([b(1),(c(1)+a(1))/2],[b(2),(c(2)+a(2))/2],...
                            'Color',[0.65 0.65 0.65],'LineWidth',0.1);
                        plot([c(1),(a(1)+b(1))/2],[c(2),(a(2)+b(2))/2],...
                            'Color',[0.65 0.65 0.65],'LineWidth',0.1);
                end
            case 4 % Square
                switch Q.Properties.Type
                    case {'Gauss-Legendre product','Gauss-Lobatto product'}
                        ax.XAxis.TickValues = unique(Q.Points(:,1));
                        ax.YAxis.TickValues = unique(Q.Points(:,2));
                    case 'nonproduct'
                        ax.XAxis.TickValues = [-1 0 1];
                        ax.YAxis.TickValues = [-1 0 1];
                        plot([-1,1],[-1,1],'Color',[0.65 0.65 0.65],'LineWidth',0.1)
                        plot([-1,1],[ 1,-1],'Color',[0.65 0.65 0.65],'LineWidth',0.1)
                end
                % Slightly expand domain to avoid points on boudary being
                % marked as outside
                domain.vertices(domain.vertices<0) = ... 
                    domain.vertices(domain.vertices<0) - eps;
                domain.vertices(domain.vertices>0) = ... 
                    domain.vertices(domain.vertices>0) + eps;             
                plot([-0.075,0.075],[0 0],'k-',[0,0],[-0.075,0.075],'k-')
                plot([-1,-0.925],[0 0],'k-',[0.925,1],[0 0],'k-',[-1,1],[0,0],'k:')
                plot([0,0],[-1,-0.925],'k-',[0,0],[0.925,1],'k-',[0,0],[-1,1],'k:')
                ax.Color = [0.831,0.816,0.784];
                ax.XAxisLocation = 'top';
                ax.TickDir = 'both';
        end
        ax.XTickLabel = {};
        ax.YTickLabel = {};
        if strcmpi(PointLabels,'on')
            text(Q.Points(:,1)+0.04, Q.Points(:,2)+0.03,...
                cellstr(num2str((1:n)')),...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','left',...
                'FontWeight','bold',...
                'Color','blue');
        end
        % Identify points with positive weights
        P = Q.Weights>0;
        % Identify points inside domain
        I = isinterior(polyshape(domain.vertices(:,1),domain.vertices(:,2)),...
            Q.Points(:,1),Q.Points(:,2));        
        % Plot points (colored according to quality)
        plot(Q.Points(P,1),Q.Points(P,2),'Marker','o',...
            'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w',...
            'LineStyle','none');
        plot(Q.Points(~P,1),Q.Points(~P,2),'Marker','o',...
            'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r',...
            'LineStyle','none');      
        plot(Q.Points(I,1),Q.Points(I,2),'Marker','o',...
            'MarkerSize',4.5,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue,...
            'LineStyle','none'); 
        plot(Q.Points(~I,1),Q.Points(~I,2),'Marker','o',...
            'MarkerSize',4.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'LineStyle','none');   
        daspect([1 1 1]);
        axis([...
            min([domain.vertices(:,1); Q.Points(:,1)]),...
            max([domain.vertices(:,1); Q.Points(:,1)]),...            
            min([domain.vertices(:,2); Q.Points(:,2)]),...
            max([domain.vertices(:,2); Q.Points(:,2)])]);
        grid off; axis off; box off;
end

end

