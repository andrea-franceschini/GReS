function ref_pts = nod2ref(pts,a,b)
% return pts coordinates in the reference space given extreme nodes

ref_pts = (2/(b-a))*(pts - (a+b)/2);

end

