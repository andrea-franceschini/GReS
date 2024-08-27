function d = computeChannelSize(p,msh1,msh2)
% Given a pair of surface mesh, compute a reasonable value for the channel
% opening
% nseg: number of segments picked for comparison
d1 = msh1.coordinates - p;
d2 = msh2.coordinates - p;

l1 = vecnorm(d1,2,2);
l2 = vecnorm(d2,2,2);

d = min(l1)+min(l2);


% compute angle between each pair of distances
% l1 = vecnorm(d1,2,2);
% l2 = vecnorm(d2,2,2);
% [~,id1] = sort(l1);
% [~,id2] = sort(l2);
% n1 = min(numel(id1),nseg);
% n2 = min(numel(id2),nseg);
% id1 = id1(1:n1);
% id2 = id2(1:n2);
% l1 = l1(id1); % pick only nseg distances
% l2 = l2(id2); % pick only nseg distances
% d1 = d1(id1,:);
% d2 = d2(id2,:);
% 
% L = l1'+l2;
% ang = zeros(n1,n2);
% % compute channel angles
% for i = 1:n1
%    for j = 1:n2
%       a = d1(i,:); 
%       b = d2(j,:);
%       ang(i,j) = acos((a*b')/(norm(a)*norm(b)));
%    end
% end
% ang = ang*180/pi;
% ind = all([ang(:)>=90 ang(:)<=200],2);
% L = L(:);
% L = L(ind);
% % keep matrix of valid channels and find smalles
% d = min(L);
end

