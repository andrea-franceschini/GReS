function  genGrain(boxSize,spheres)
% this function generate a portion of a gmsh file with spheroidal grains
% located inside a box of specified size
fID = fopen('grains.dat','w');
% convert array in comma separated string
bstr = makeString(boxSize);
fprintf(fID,'Box(1) = {%s};\n',bstr);
for i = 1:size(spheres,1)
   s = makeString(spheres(i,:));
   fprintf(fID,'Sphere(%i) = {%s};\n',i+1,s);
end

k = size(spheres,1)+2; % volume counter
% boolean differences
for i = 1:size(spheres,1)-1
   fprintf(fID,'BooleanIntersection(%i) = { Volume{1};}{ Volume{%i}; Delete;};\n',k,i+1);
   k = k+1;
end
fprintf(fID,'BooleanIntersection(%i) = { Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',k,i+2);
