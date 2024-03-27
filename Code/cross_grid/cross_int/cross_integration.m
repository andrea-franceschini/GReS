%%
%TESTING CROSS INTEGRATION OVER A PATCH OF NON MATCHING 1D ELEMENTS
close all; clear

% master surfaces node position along 1D axis
master = (linspace(-3,2,5))';
slave = (linspace(-3,2,14))';
% matrix in which to store results
matrix = zeros(length(master), length(slave));

nInt = 5; % number of RBF interpolation points for each element (min. 2)
nGP = 6; % number of Gauss points (value admitted by Gauss class) 

% build a topology based on nodes position
% assume node in ordered sequence
mastertop = build_topol(master);
slavetop = build_topol(slave);

% build connectivity between elements - check for intersection
elem_connectivity = zeros(size(mastertop,1), size(slavetop,1));

for i = 1:size(mastertop,1)
    a = master(mastertop(i,1));
    b = master(mastertop(i,2));
    % loop trough master element to find connectivity
    for j = 1:size(slavetop,1)
        c = slave(slavetop(j,1));
        d = slave(slavetop(j,2));
        if ~any([a>d,c>b])
            % intersecting
            elem_connectivity(i,j) = 1;
        end
    end
end


for i = 1:length(master)
    % get element of the support
    master_elems = find(mastertop == nodeID);
    % get shape function values and interpolation points coordinate
    [fInt, ptsInt] = evalSF(i, master_elems, nInt, mastertop, master);
    % compute RBF weight of interpolant
    for ii = 1:length(ptsInt)
        r = abs(max(ptsInt) - min(ptsInt));
        d = ptsInt - ptsInt(ii);
        d = sqrt(d.^2);
        rbf = pos(1-d./r).^4.*(1+d./r);
        % continue but using full matrix, not sparse
    end
    % copy code from local pc test

    % get elements on slave surface connected to master support of node i
    slave_elems = elem_connectivity(master_elems,:);
    % from list of slave elems, compute interpolant on GaussPoints (in real
    % space)
    
    % loop trough slave_elems, and popoulate matrix entries




end
