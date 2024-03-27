function real_pts = ref2nod(pts,a,b)
% return pts coordinates in the real space given extreme nodes
% a = [x1,y1]    b = [x2,y2];
real_pts(:,1) = ((b(1)-a(1))/2)*pts + (a(1)+b(1))/2;
real_pts(:,2) = zeros(length(pts),1);
end

