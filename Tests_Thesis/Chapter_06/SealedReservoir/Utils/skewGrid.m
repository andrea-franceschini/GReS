function skewGrid(grid, skew)
% SKEWGRID Apply simultaneous planar skews to a structured box grid.
%
%   skewGrid(grid, skew)
%
%   skew is a struct array with fields:
%     - surfTag : boundary surface tag, from 1 to 6
%     - normal  : requested outward normal of the skewed face
%
%   Surface tags:
%     1 bottom zmin
%     2 top    zmax
%     3 south  ymin
%     4 north  ymax
%     5 west   xmin
%     6 east   xmax

if ~grid.isStructured()
    error("skewGrid:NonStructuredGrid", ...
      "skewGrid is only valid for structured grids.");
end

X0 = grid.coordinates;

tol = 1e-10 * max(1, norm(X0, "fro"));

if ~isAxisAlignedBoxGrid(X0, tol)
    error("skewGrid:AlreadySkewed", ...
      "Grid is not an axis-aligned structured box. Apply all skews simultaneously from the original structured grid.");
end

if ~isstruct(skew) || ~all(isfield(skew, ["surfTag", "normal"]))
    error("skewGrid:InvalidInput", ...
      "Input skew must be a struct array with fields surfTag and normal.");
end

surfTags = [skew.surfTag];

if numel(unique(surfTags)) ~= numel(surfTags)
    error("skewGrid:RepeatedSurface", ...
      "Each surface can be skewed only once in the same call.");
end

if any(surfTags < 1 | surfTags > 6 | abs(surfTags - round(surfTags)) > 0)
    error("skewGrid:InvalidSurface", ...
      "Surface tags must be integers from 1 to 6.");
end

xmin = min(X0(:, 1));
xmax = max(X0(:, 1));
ymin = min(X0(:, 2));
ymax = max(X0(:, 2));
zmin = min(X0(:, 3));
zmax = max(X0(:, 3));

if any([xmax - xmin, ymax - ymin, zmax - zmin] <= 0)
    error("skewGrid:DegenerateBox", ...
      "The structured grid has degenerate box dimensions.");
end

dX = zeros(size(X0));

for k = 1:numel(skew)

    tag = skew(k).surfTag;
    n = skew(k).normal(:).';

    if numel(n) ~= 3 || ~isreal(n) || any(~isfinite(n)) || norm(n) == 0
        error("skewGrid:InvalidNormal", ...
          "Normal must be a finite real non-zero 3-vector.");
    end

    n = n / norm(n);

    [e, p0, w, Xface] = localFaceData(tag, X0, ...
      xmin, xmax, ymin, ymax, zmin, zmax);

    c = dot(e, n);

    if c <= 0
        error("skewGrid:IncompatibleNormal", ...
          "Requested normal for surface %d must have positive projection on the original outward normal.", tag);
    end

    if c < 0.25
        error("skewGrid:ExcessiveSkew", ...
          "Requested normal for surface %d is too close to tangential to the original face.", tag);
    end

    % Project the corresponding boundary-face point onto the requested plane.
    % For example, for the east face this uses Xface = [xmax, y, z],
    % not the actual interior point X = [x, y, z].
    alpha = ((p0 - Xface) * n.') / c;

    dX = dX + w .* (alpha .* e);

end

grid.coordinates = X0 + dX;

end

function tf = isAxisAlignedBoxGrid(X, tol)

x = uniquetol(X(:, 1), tol);
y = uniquetol(X(:, 2), tol);
z = uniquetol(X(:, 3), tol);

tf = numel(x) * numel(y) * numel(z) == size(X, 1);

end

function [e, p0, w, Xface] = localFaceData(tag, X, ...
  xmin, xmax, ymin, ymax, zmin, zmax)

Xface = X;

switch tag
    case 1 % bottom, zmin
        e  = [0 0 -1];
        p0 = [mean([xmin xmax]), mean([ymin ymax]), zmin];
        s  = (zmax - X(:, 3)) / (zmax - zmin);
        Xface(:, 3) = zmin;

    case 2 % top, zmax
        e  = [0 0 1];
        p0 = [mean([xmin xmax]), mean([ymin ymax]), zmax];
        s  = (X(:, 3) - zmin) / (zmax - zmin);
        Xface(:, 3) = zmax;

    case 3 % south, ymin
        e  = [0 -1 0];
        p0 = [mean([xmin xmax]), ymin, mean([zmin zmax])];
        s  = (ymax - X(:, 2)) / (ymax - ymin);
        Xface(:, 2) = ymin;

    case 4 % north, ymax
        e  = [0 1 0];
        p0 = [mean([xmin xmax]), ymax, mean([zmin zmax])];
        s  = (X(:, 2) - ymin) / (ymax - ymin);
        Xface(:, 2) = ymax;

    case 5 % west, xmin
        e  = [-1 0 0];
        p0 = [xmin, mean([ymin ymax]), mean([zmin zmax])];
        s  = (xmax - X(:, 1)) / (xmax - xmin);
        Xface(:, 1) = xmin;

    case 6 % east, xmax
        e  = [1 0 0];
        p0 = [xmax, mean([ymin ymax]), mean([zmin zmax])];
        s  = (X(:, 1) - xmin) / (xmax - xmin);
        Xface(:, 1) = xmax;
end

s = max(0, min(1, s));

% Linear blending is safer than smoothstep here.
% It preserves non-zero layer thickness near the skewed face.
w = s;

end
