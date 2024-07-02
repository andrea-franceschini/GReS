function  genGrain(fG,fW,boxSize,spheres)
% this function generate a portion of a gmsh file with spheroidal/ellipsoidal grains
% located inside a box of specified size
fIDG = fopen(fG,'w');
fIDW = fopen(fW,'w');
% convert array in comma separated string
bstr = makeString(boxSize);
fprintf(fIDG,'Box(1) = {%s};\n',bstr);
fprintf(fIDW,'Box(1) = {%s};\n',bstr);
for i = 1:numel(spheres)
   if numel(spheres{i})<5
      s = makeString(spheres{i});
      fprintf(fIDG,'Sphere(%i) = {%s};\n',i+1,s);
      fprintf(fIDW,'Sphere(%i) = {%s};\n',i+1,s);
   else
      s = makeString(spheres{i}(1:4));
      fprintf(fIDG,'Sphere(%i) = {%s};\n',i+1,s);
      fprintf(fIDW,'Sphere(%i) = {%s};\n',i+1,s);
      d = makeString(spheres{i}(5:end));
      fprintf(fIDG,'Dilate{%s} = {Volume{%i};};\n',d,i+1);
      fprintf(fIDW,'Dilate{%s} = {Volume{%i};};\n',d,i+1);
      %Dilate {2, 1, 1} { Volume{3};}
   end
end

k = numel(spheres)+2; % volume counter
% boolean interesection for grain file
for i = 1:numel(spheres)-1
   fprintf(fIDG,'BooleanIntersection(%i) = { Volume{1};}{ Volume{%i}; Delete;};\n',k,i+1);
   fprintf(fIDW,'BooleanDifference{ Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',i+1);
   k = k+1;
end
fprintf(fIDG,'BooleanIntersection(%i) = { Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',k,i+2);
fprintf(fIDW,'BooleanDifference{ Volume{1}; Delete;}{ Volume{%i}; Delete;};\n',i+2);
end