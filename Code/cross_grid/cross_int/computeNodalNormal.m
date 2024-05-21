function n_n = computeNodalNormal(topol,coord)
% Return vector n of weighted normals

% compute tangent of each element
t = [coord(topol(:,2),1) - coord(topol(:,1),1), coord(topol(:,2),2) - coord(topol(:,1),2)];
t = t./sqrt((t(:,1).^2 + t(:,2).^2));
n = [t(:,2), -t(:,1)]; % chosen arbitrarily, later this will be checked

% compute length of each element
l = vecnorm([coord(topol(:,1),1) - coord(topol(:,2),1), coord(topol(:,1),2) - coord(topol(:,2),2)],2,2);
n_n = zeros(length(coord),2);

for i = 1:length(coord)
    el = find(any((ismember(topol,i)),2));
    if length(el) > 1
        n_n(i,:) = (l(el(1))*n(el(1),:) + l(el(2))*n(el(2),:))/(l(el(1))+l(el(2)));
    elseif length(el) == 1
        n_n(i,:) = n(el,:);
    elseif isempty(el)
        n_n(i,:) = 0;
    end
end


