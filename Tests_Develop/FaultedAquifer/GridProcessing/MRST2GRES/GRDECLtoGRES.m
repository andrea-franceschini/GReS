function [msh, conn] = GRDECLtoGRES(grdecl,G,rat)

% build a GRES mesh object starting from an eclipse grid
% leverages MRST routines

% mesh - GReS Mesh() object
% surfTag 1 -> slave faces
% surfTag 2 -> master faces

% G: processed eclipse grid

% rat: minimum ratio between vertical edges of the cell and x-y extent
% (avoid to thin cells)

tol = rat*mean(mean(G.cells.cpgeometry.extent(:,[1 2]),1));

% rerun MRST processor with tolerance to merge points on a pillar
% thid modified version return also the cell connecitivity along the faulted
% pillars in i and j direction
G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',tol);

% use face tag and neighbors to build topology of mortar faces
i = G.faces.tag == 1; % faulted faces

is_i_oriented = i & G.faces.cellTags(:,1) == 1;

is_j_oriented = i & G.faces.cellTags(:,1) == 3;

% distinguish neighbor pairs of master and slave side
n_i = G.faces.neighbors(is_i_oriented,:);
n_j = G.faces.neighbors(is_j_oriented,:);

% remove boundary faces
n_i = n_i(all(n_i,2),:);
n_j = n_j(all(n_j,2),:);



% according to MRST face orientation and mortar pattern:
% n_i(:,1) -> slave n_i(:,2) -> master
% n_j(:,1) -> master n_i(:,2) -> slave
% this should avoid cross point issues


% get face topology using preserved corner point nodes
%                            +-----------> I
%                           /|
%                          / |     1 --------- 2
%                         /  |    /|          /|   Column numbers mapped
%                      J v   |   / |         / |   to vertices in ECLIPSE's
%                            |  3 --------- 4  |   default right-handed
%                          K v  |  |        |  |   coordinate system
%                               |  5 -------|- 6   (origin is top, left,
%                               | /         | /    back vertex with Z-axis
%                               |/          |/     pointing down.)
%                               7 --------- 8

% slave i : 2-4-6-8
% master i: 1-3-5-7
% slave j = 1-2-5-6
% master j = 4-3-8-7

% get face topology
[s_i,~,s_ii] = unique(n_i(:,1)); % cell of slave side in i direction
[m_i,~,m_ii] = unique(n_i(:,2)); % cell of master side in i direction
[s_j,~,s_jj] = unique(n_j(:,2)); % cell of slave side in j direction
[m_j,~,m_jj] = unique(n_j(:,1)); % cell of master side in j direction

[nsi,nmi,nsj,nmj]  = deal(numel(s_i),numel(m_i),numel(s_j),numel(m_j));


nS = sum([nsi nsj nmi nmj]);

cells = G.cells.cpnodes;

% build face topology matrix
surf_si = cells(s_i,[2 4 6 8]);
surf_mi = cells(m_i,[3 1 7 5]);
surf_sj = cells(s_j,[1 2 5 6]);
surf_mj = cells(m_j,[4 3 8 7]);

m_ii = m_ii + nsi;
s_jj = s_jj + nsi + nmi;
m_jj = m_jj + nsi + nmi + nsj;

% reconstruct connectivity matrix using indices of surfaces
conn = sparse([m_ii;m_jj],[s_ii;s_jj],1,nS,nS);
conn(conn~=0) = true;

% assign surfaceTags to distinguish master and slave faces
% 1 -> slave
% 2 -> master


% duplicate coinciding nodes in faulted faces

% if two nodes coincide, the node of the master receives a new node index.
% the cell indices node must be modified accordingly
[x_i,~,id_i] = intersect(surf_si,surf_mi);
[x_j,~,id_j] = intersect(surf_sj,surf_mj);

% update new node numbers in intersected nodes of i master faces
new = (1:numel(x_i))';
surf_mi(id_i) = new + G.nodes.num;

% update new node numbers in intersected nodes of j master faces
% caution: do not duplicate nodes that have already been duplicated on i
% faces

% find entries of x_j that belong to x_i
is_already_duplicated = ismembc(x_j,x_i);
% remove corresponding logical index from id_j
id_j(is_already_duplicated) = [];
new = (1:sum(~is_already_duplicated))';

surf_mj(id_j) = new + G.nodes.num + numel(x_i);


% reupdate cell topology
cells(m_i,[3 1 7 5]) = surf_mi;
cells(m_j,[4 3 8 7]) = surf_mj;

surfaces = [surf_si; ...
  surf_sj; ...
  surf_mi;...
  surf_mj];

surfaceTag = [ones(nsi+nsj,1);...
  2*ones(nmi+nmj,1)];

% define GReS object
msh = Mesh();

msh.nDim = 3;
nN = G.nodes.num + numel(x_i) + numel(x_j(~is_already_duplicated));
msh.nNodes = nN;

% split collapsed hexahedron into valid tetrahedron
cells = cells(:,[1 2 4 3 5 6 8 7]);         % valid VTK ordering
cells = hex2tet(cells);

% restore vtk ordering
msh.cells = cells;
msh.nCells = size(msh.cells,1);
msh.cellTag = ones(msh.nCells,1);
msh.cellNumVerts = sum(msh.cells > 0,2);
msh.cellVTKType = zeros(msh.nCells,1);
msh.cellVTKType(msh.cellNumVerts == 4) = 10;
msh.cellVTKType(msh.cellNumVerts == 8) = 12;
% assume for now mesh made of only hexahedra (even if nodes are duplicated)


surfaces = surfaces(:,[1 2 4 3]);   % use vtk indexing
% transform degenerate faces in triangles

% spot degenerate faces
[ss,is] = sort(surfaces,2);
diff_surf = diff(ss,1,2);
is_surf_degenerate = ~all(diff_surf,2);
[dup_vert,~] = find(diff_surf(is_surf_degenerate,:)'==0);
id = find(is_surf_degenerate);
k = 0;
for i = id'
  % remove index and make 0 the last element 
  k = k+1;
  p = [1 2 3 4] == is(i,dup_vert(k)+1);
  surfaces(i,:) = surfaces(i,[find(~p) find(p)]);
  surfaces(i,4) = 0;
end



msh.surfaces = surfaces;
msh.surfaceTag = surfaceTag;
msh.nSurfaces = size(msh.surfaces,1);
msh.surfaceNumVerts = sum(msh.surfaces>0,2);
msh.surfaceVTKType = zeros(msh.nSurfaces,1);
msh.surfaceVTKType(msh.surfaceNumVerts == 3) = 5;
msh.surfaceVTKType(msh.surfaceNumVerts == 4) = 9;



% update coordinates including new nodes
msh.coordinates = G.nodes.coords;
msh.coordinates = [msh.coordinates;...
                    msh.coordinates(x_i,:);...
                    msh.coordinates(x_j(~is_already_duplicated),:)];

% better safe then sorry
% remove left over nodes from the coordinates list

nodeList = unique(msh.cells);

assert(numel(nodeList)==max(msh.cells,[],"all"),"Node renumbering required");
msh.coordinates = msh.coordinates(nodeList,:);
msh.nNodes = size(msh.coordinates,1);
  

end




