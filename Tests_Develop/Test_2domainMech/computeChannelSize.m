function d = computeChannelSize(p,msh1,msh2)
% Given a pair of surface mesh, compute a reasonable value for the channel
% opening
nseg = 20; % number of segments picked for comparison
d1 = msh1.coordinates - p;
d2 = msh2.coordinates - p;
% compute angle between each pair of distances
l1 = vecnorm(d1,2,2);
l2 = vecnorm(d2,2,2);
[~,id1] = sort(l1);
[~,id2] = sort(l2);
id1 = id1(1:nseg);
id2 = id2(1:nseg);
l1 = l1(id1); % pick only 20 distances
l2 = l2(id2); % pick only 20 distances
d1 = d1(id1,:);
d2 = d2(id2,:);

L = l1'+l2;
ang = zeros(nseg);
% compute channel angles
for i = 1:nseg
   for j = 1:nseg
      a = d1(i,:); 
      b = d2(j,:);
      ang(i,j) = acos((a*b')/(norm(a)*norm(b)));
   end
end
ang = ang*180/pi;
ind = all([ang(:)>=160 ang(:)<=200],2);
L = L(:);
L = L(ind);
% keep matrix of valid channels and find smalles
d = min(L);
end

