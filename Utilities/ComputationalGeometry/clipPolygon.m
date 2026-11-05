function polyClip = clipPolygon(poly1,poly2)

% Clip 2D polygons
% poly1,poly2: N x 2 coordinates of vertices
% Uses Sutherland-Hodgman clipping algorithm
% Algortihm assumes polygons are 2D - convex - CCW ordered
% Return clipped polygon

% order points in CCW order
p = mxOrderPointsCCW(poly1);
poly1 = poly1(p,:);

p = mxOrderPointsCCW(poly2);
poly2 = poly2(p,:);

% clip using mex function
polyClip = mxPolygonClip(poly1,poly2);

end
