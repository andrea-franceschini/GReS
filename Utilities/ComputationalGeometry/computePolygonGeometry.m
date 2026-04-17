function varargout = computePolygonGeometry(varargin)
%COMPUTEPOLYGONGEOMETRY Compute area, centroid, and normal of 2D/3D polygons.
%
%  [A,C,N] = COMPUTEPOLYGONGEOMETRY(P)
%  [A,C,N] = COMPUTEPOLYGONGEOMETRY(P,NORMAL)
%  X       = COMPUTEPOLYGONGEOMETRY(...,MODE)
%
%  Thin MATLAB wrapper around the mex function MXPOLYGONGEOMETRY.
%
%  Local syntax
%  ------------
%  P is an N-by-2 (or N-by-3) array containing the vertices of one polygon.
%
%    [A,C,N] = COMPUTEPOLYGONGEOMETRY(P)
%      returns area A, centroid C, and, in 3D, unit normal N.
%
%    [A,C,N] = COMPUTEPOLYGONGEOMETRY(P,NORMAL)
%      same as above, but enforces the orientation using the supplied
%      3-vector NORMAL. This is useful in 3D when the polygon orientation
%      must match an external convention.
%
%    X = COMPUTEPOLYGONGEOMETRY(...,MODE)
%      returns only the quantity selected by MODE:
%        "area"      -> scalar area
%        "centroid"  -> 1-by-dim centroid
%        "normal"    -> 1-by-3 normal (3D only)
%        "geometry"  -> same as requesting [A,C,N]
%
%  Batch syntax
%  ------------
%  P may also contain multiple polygons concatenated row-wise, with NVERT
%  specifying the number of vertices of each polygon.
%
%    [A,C,N] = COMPUTEPOLYGONGEOMETRY(Pflat,NVERT)
%    [A,C,N] = COMPUTEPOLYGONGEOMETRY(Pflat,NVERT,NORMALS)
%    X       = COMPUTEPOLYGONGEOMETRY(...,MODE)
%
%  where:
%    Pflat    : sum(NVERT)-by-dim array of stacked polygon vertices
%    NVERT    : vector with one entry per polygon
%    NORMALS  : nPoly-by-3 array of prescribed normals (3D only)
%
%  In batch mode, outputs are returned row-wise, one row per polygon:
%    A : nPoly-by-1
%    C : nPoly-by-dim
%    N : nPoly-by-3  (3D only)
%
%  Notes
%  -----
%  If normals are not supplied, polygons are assumed to already be ordered
%  consistently along the perimeter, with no self-crossing edges.
%
%  If a normal is supplied for a 3D polygon, the mex routine uses it to
%  enforce a consistent orientation; this is useful when points are
%  coplanar but their ordering may not follow the desired convention.
%
%  See also MXPOLYGONGEOMETRY.


[areas,centers,normals] = mxPolygonGeometry(varargin{:});


varargout{1} = areas;
varargout{2} = centers;

if nargout > 2
  varargout{3} = normals;
end

end

