function cells = hex2tet(cells)


% test nodes 5 6 7 8 and 1 2 3 4 (check degenerate cases based on indices)
deg_nodes = cells(:,[5 6 7 8]) ==  cells(:,[1 2 3 4]); 

deg_cells = all(deg_nodes,2);

assert(~any(deg_cells),...
  'detected %i pinched cells\n',sum(deg_cells))

% find cells according to degenerate cases to process them individually
df1 = reshape(find(sum(deg_nodes,2)==3),1,[]);
df2 = reshape(find(sum(deg_nodes,2)==2),1,[]);
df3 = reshape(find(sum(deg_nodes,2)==1),1,[]);

new_tet = zeros(0,4);

% 3 pinched node
for i = df1
  k = find(~deg_nodes(i,:));
  nodes = cells(i,:);
  new_tet = [new_tet;...
    nodes([k v(k+2) v(k+3) k+4]);...
    nodes([k v(k+1) v(k+2)  k+4])];
end

% 2 pinched nodes
for i = df2
  nodes = cells(i,:);
  k = find(~deg_nodes(i,:));
  if abs(k(2)-k(1))==2
    k1 = k(1); k2 = k(2);
    % opposite pinchde nodes
    new_tet = [new_tet;...
      nodes([k1 v(k1+1) v(k2+1) k1+4]);...
      nodes([v(k1+1) k2 v(k2+1) k2+4])];
  else
    % fixed combination: the smaller k and its opposite are always
    % present
    if k(1)==1 && k(2)==4
      k = fliplr(k);
    end

    k1 = k(1); k2 = k(2);
    nList = [v(k1-1) k1 v(k1+2) k1+4;...
             v(k2-1) k2 v(k2+1) k2+4;...
             k1 v(k1+2) k2+4 k1+4];
    assert(numel(unique(nList))==6,['Wrong numbering' ...
      ' of tetrahedron subdivision for hexa %i\n'],i)
    new_tet = [new_tet;
               nodes(nList)];
  end
end

% 1 pinched node
% 4 tetra with all diagonals built from the pinched node
for i = df3
  nodes = cells(i,:);
  k = find(deg_nodes(i,:));
  new_tet = [new_tet;...
    nodes([k v(k+1) v(k+2) v(k+1)+4]);...
    nodes([k v(k+1)+4 v(k+2) v(k+2)+4]);...
    nodes([k v(k+2) v(k+3) v(k+3)+4]);...
    nodes([k v(k+2) v(k+3)+4 v(k+2)+4])];
end

% remove degenerated cells and appen new tetra
cells([df1';df2';df3'],:) = [];
new_tet = [new_tet zeros(size(new_tet,1),4)];
cells = [cells; new_tet];

end

function b = v(a)
 % sum node index in loop [1 2 3 4]
 v = [1 2 3 4];
 b = v(mod(a-1, 4) + 1);
end

