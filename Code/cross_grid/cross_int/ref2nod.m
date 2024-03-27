function real_pts = ref2nod(pts,a,b)
% return pts coordinates in the real space given extreme nodes

real_pts = ((b-a)/2)*pts + (a+b)/2;

end

