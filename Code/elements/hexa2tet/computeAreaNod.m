function aNod = computeAreaNod(mesh,face)
%COMPUTEAREANOD Summary of this function goes here
%   Detailed explanation goes here
aNod = zeros(mesh.nNodes,1);
for el=1:mesh.nSurfaces
    % Get the right material stiffness for each element
    switch mesh.surfaceVTKType(el)
        case 5 % Triangle
           top = mesh.surfaces(el,:);
           aNod(top) = aNod(top) + face.computeAreaTri(el)/3;
            %             end
        case 9 % Quadrilateral
           top = mesh.surfaces(el,:);
           aNod(top) = aNod(top) + face.computeAreaQuad(el)/4;
    end
end

