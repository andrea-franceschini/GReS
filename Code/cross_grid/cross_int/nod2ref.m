function ref_pts = nod2ref(pts,a,b)
% return pts coordinates in the reference space given extreme nodes
ref_pts = -1 + 2*(pts(:,1)-a(1))/(b(1)-a(1));
end


% just linear interpolation

