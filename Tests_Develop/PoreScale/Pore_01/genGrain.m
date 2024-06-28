function  genGrain(fG,fW,boxSize,spheres)
% this function generate a portion of a gmsh file with spheroidal grains
% located inside a box of specified size
fIDG = fopen(fG,'w');
fIDW = fopen(fW,'w');
% convert array in comma separated string
bstr = makeString(boxSize);
fprintf(fIDG,'Box(1) = {%s};\n',bstr);
fprintf(fIDW,'Box(1) = {%s};\n',bstr);
for i = 1:size(spheres,1)
   s = makeString(spheres(i,:));
   fprintf(fIDG,'Sphere(%i) = {%s};\n',i+1,s);
   fprintf(fIDW,'Sphere(%i) = {%s};\n',i+1,s);
end

k = size(spheres,1)+2; % volume counter
% boolean interesection for grain file
for i = 1:size(spheres,1)-1
   fprintf(fIDG,'BooleanIntersection(%i) = { Volume{1};}{ Volume{%i}; Delete;};\n',k,i+1);
   fprintf(fIDW,'BooleanDifference{ Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',i+1);
   k = k+1;
end
fprintf(fIDG,'BooleanIntersection(%i) = { Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',k,i+2);
fprintf(fIDW,'BooleanDifference{ Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',i+2);
end