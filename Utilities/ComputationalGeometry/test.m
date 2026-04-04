poly1 = [0 0 0; 1 0 0; 0 1 0];
poly2 = [0 0 1; 1 0 1; 1 1 1; 0 1 1];

P = [poly1; poly2];
nVert = [3; 4];

[A,c,n] = computePolygonGeometry(P,nVert);

poly1 = poly1(:,1:end-1);
poly2 = poly2(:,1:end-1);
clip = clipPolygon(poly1,poly2);