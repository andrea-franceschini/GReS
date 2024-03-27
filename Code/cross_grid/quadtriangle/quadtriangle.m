function Q = quadtriangle(d,varargin)

%% QUADTRIANGLE Quadrature weights and points for a triangle
%    Q = QUADTRIANGLE(d) returns the weights and points of a Gauss-Jacobi 
%    product rule of degree d for numerically integrating over a (default) 
%    triangle defined by the vertices (-1,-1),(1,-1) and (-1,1). The output 
%    is stored as a structure Q with the following fields:
%       »  Q.Points = the n quadrature points in an array of size n×2, 
%          where Q.Points(:,1) and Q.Points(:,2) are the x and y 
%          coordinates, respectively.
%       »  Q.Weights = the corresponding quadrature weights stored in a 
%          column vector of length n.
%       »  Q.Properties = contains (sub)fields describing the .Degree (the 
%          input value d), the .Type (see description below) and the 
%          .Domain of integration of the quadrature rule.
%  
%    Q = QUADTRIANGLE(d,Name,Value) same as above but with additional 
%    control over the type of rules returned. Options include specification 
%    of the following Name-Value pair arguments: 
%    ______________________________________________________________________
%    'Domain' -- Triangular domain
%    [ -1 -1; 1 -1; -1 1 ] (default) | [ x1 y2; x2 y2; x3 y3 ] | [ ] 
%    ----------------------------------------------------------------------     
%    The triangular domain over which the weights and points are defined,
%    specified as a 3×2 array of the (x,y) vertices of the triangle, i.e.,
%    T = [x1 y2; x2 y2; x3 y3], or as an empty array. 
% 
%    In the case of the latter option, the returned quadrature points 
%    consist of the first two barycentric (or area) coordinates, and the 
%    returned quadrature weights are scaled such that they sum to 1. (Note: 
%    The third barycentric coordinate for the points can easily be computed 
%    as 1-Q.Points(:,1)-Q.Points(:,2).)
%    ______________________________________________________________________
%    'Type' -- Type of quadrature rule
%    'product' (default) | 'nonproduct'
%    ----------------------------------------------------------------------  
%    The type of quadrature rule may be specified as 'product' (the
%    default), in which case the points and weights of a Gauss-Jacobi
%    product rule are returned, or as 'nonproduct'. Product rules can be of
%    arbitrary degree, while nonproduct rules are currently available up to
%    degree 25.
%
%    The following additional Name-Value pairs are available for nonproduct 
%    rules:
%    ______________________________________________________________________
%    'Symmetry' -- Symmetry of the quadrature rule
%    'full' (default) | 'allowAsymmetric'
%    ----------------------------------------------------------------------  
%    The 'full' (default) value requires the quadrature rule to be fully 
%    symmetric. The value 'allowAsymmetric' allows asymmetric quadrature 
%    rules of lower point count than the full symmetry rules to be returned 
%    when available. (Note: Some quadrature rules possess so-called partial
%    symmetry, though no distinction is made here between asymmetric and
%    partially symmetric rules, i.e., a rule is classified as either fully
%    symmetric or asymmetric.)
%    ______________________________________________________________________
%    'Weights' -- Quadrature weight conditions
%    'positive' (default) | 'allowNegative'
%    ---------------------------------------------------------------------- 
%    The 'postive' (default) value requires all quadrature weights to be 
%    positive. The value 'allowNegative' allows quadrature rules with nega-
%    tive weight(s) of lower point count than the positive weight rules to
%    be returned when available.
%    ______________________________________________________________________
%    'Points' -- Location of quadrature points
%    'inside' (default) | 'allowOutside'
%    ---------------------------------------------------------------------- 
%    The 'inside' (default) value requires all points to be inside the tri-
%    angle. The value 'allowOutside' allows quadrature rules with point(s)  
%    outside the triangle of lower point count than the quadrature rules 
%    with all points inside to be returned when available.
%    ______________________________________________________________________                 
%
%    Note: The default Name-Value pairs for nonproduct rules are set such 
%    that the minimal-point, fully symmetric so-called PI rule (that is,  
%    all weights are Positive and all points are Inside the triangle) of 
%    degree d is returned. Fully symmetric PI rules are generally preferred 
%    among practitioners.
%
%    See also quadsquare, quadGaussJacobi, quadGaussLobatto
%

%% Validate data passed to function

vd = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0});
vD = @(x)validateattributes(x,{'numeric'},{}); 
vT = @(x)validateattributes(x,{'char'},{'nonempty'});
ip = inputParser;
ip.addRequired('d',vd);
ip.addParameter('Domain',[ -1,-1; 1,-1; -1, 1 ],vD);
ip.addParameter('Type','product',vT);
ip.addParameter('Symmetry','full')
ip.addParameter('Weights','positive'),
ip.addParameter('Points','inside');
ip.parse(d,varargin{:}); ip.Results; 
Domain = ip.Results.Domain; Type = ip.Results.Type; 
Symmetry = ip.Results.Symmetry; Weights  = ip.Results.Weights; Points = ip.Results.Points; 

%% Compute quadrature weights and points
switch Type
    case {'product','Product'}
        %% Product Rules
        Qn = quadTriangleProduct(d);
    case {'nonproduct','Nonproduct'}
        %% Nonproduct rules
        if ~any(strcmpi(Symmetry,{'full','allowAsymmetric'}))
            error(['The value of ''Symmetry'' is invalid. Expected input to be:',...
                ' ''full'' (default) or ''allowAsymmetric''.']);
        end
        if ~any(strcmpi(Weights,{'positive','allowNegative'}))
            error(['The value of ''Weights'' is invalid. Expected input to be:',...
                ' ''positive'' (default) or ''allowNegative''.']);
        end
        if ~any(strcmpi(Points,{'inside','allowOutside'}))
            error(['The value of ''Points'' is invalid. Expected input to be:',...
                ' ''inside'' (default) or ''allowOutside''.']);
        end
        Qn = quadTriangleNonproduct(d,Symmetry,Weights,Points);
    otherwise
        error(['The value of ''Type'' is invalid. Expected input to be:',...
            ' ''product'' or ''nonproduct''.']);
end
Q = quadTrianglePointsAndWeights(Qn,Domain);
Q.Properties.Degree = d;
Q.Properties.Type   = Type;
Q.Properties.Domain = Domain;
end

%% Supporting function(s) =================================================