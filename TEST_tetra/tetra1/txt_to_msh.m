close all;
clear;
clc;

%--------------------- FROM TXT TO MSH ------------------
fileName = 'tetra1.txt';
%--------------------------------------------------------

fid = fopen(fileName,'r');

while (~feof(fid))
     line = fgetl(fid);
     parts = strsplit(line,' ');
  if length(parts) > 1    
     line = parts{2};
  end
     while (isequal(line, 'Node'))
         nodeID = str2double(parts{3});
         coord_x = str2double(parts{4});
         coord_y = str2double(parts{5});
         coord_z = str2double(parts{6});
         Coords(nodeID,:) = [coord_x coord_y coord_z];
         line = fgetl(fid);
         parts = strsplit(line,' ');
         if length(parts) > 1 
         line = parts{2};
         end
     end
   
     while (isequal(line, 'Tri3'))
         surfID = str2double(parts{3});
         surf_PhyName = str2double(parts{4});
         surfTag = str2double(parts{5});
         n1 = str2double(parts{6});
         n2 = str2double(parts{7});
         n3 = str2double(parts{8});
         surfTopol(surfID,:) = [n1 n2 n3];
         line = fgetl(fid);
         parts = strsplit(line,' ');
         if length(parts) > 1 
         line = parts{2};
         end
     end
     while (isequal(line, 'Quad4'))
         surfID = str2double(parts{3});
         surf_PhyName = str2double(parts{4});
         surfTag = str2double(parts{5});
         n1 = str2double(parts{6});
         n2 = str2double(parts{7});
         n3 = str2double(parts{8});
         n4 = str2double(parts{9});
         surfTopol(surfID,:) = [n1 n2 n3 n4];
         line = fgetl(fid);
         parts = strsplit(line,' ');
         if length(parts) > 1 
         line = parts{2};
         end
     end
     while (isequal(line, 'Tetra4'))
         elemID = str2double(parts{3});
         el_PhyName = str2double(parts{4});
         cellTag = str2double(parts{5});
         n1 = str2double(parts{6});
         n2 = str2double(parts{7});
         n3 = str2double(parts{8});
         n4 = str2double(parts{9});
         elemTopol(elemID,:) = [n1 n2 n3 n4];
         line = fgetl(fid);
         parts = strsplit(line,' ');
         if length(parts) > 1 
         line = parts{2};
         end
     end
     i = 1;
     while (isequal(line, 'Hexa8'))
         elemID = str2double(parts{3});
         el_PhyName(i) = str2double(parts{4});
         cellTag(i) = str2double(parts{5});
         n1 = str2double(parts{6});
         n2 = str2double(parts{7});
         n3 = str2double(parts{8});
         n4 = str2double(parts{9});
         n5 = str2double(parts{10});
         n6 = str2double(parts{11});
         n7 = str2double(parts{12});
         n8 = str2double(parts{13});
         elemTopol(elemID,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
         line = fgetl(fid);
         parts = strsplit(line,' ');
         if length(parts) > 1 
          line = parts{2};
          i = i+1;
         end
     end

end

 fclose(fid);

%--------------------------- FILE MSH ---------------------------------
 
 [nRow,~] = size(Coords);
 nTotNodes = nRow;
 [nRow,nCol] = size(elemTopol);
 nElemNodes = nCol;
 nTotElems = nRow;
%  [nRow,nCol] = size(surfTopol);
%  nTotSurfs = nRow;
 
 if nElemNodes == 4
   codeSup = 2;
   codeEl = 4;
 elseif nElemNodes == 8
   codeSup = 3;
   codeEl = 5;
 end
 
 num = 1 : nTotNodes;
 nTotSurfs = 0;
 nTotEn = nTotElems + nTotSurfs;
%  surf_flag = linspace(2,2,nTotSurfs);
 elem_flag = linspace(2,2,nTotElems);
%  nPhyNames = length(el_PhyName)+length(surf_PhyName);
 el_PhyName = 1;
 nPhyNames = max(el_PhyName);

 
 fid = fopen('mesh.msh','w');
 fprintf(fid,'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
 fprintf(fid,'$PhysicalNames\n');
 fprintf(fid,'%g\n',nPhyNames);
 fprintf(fid,'3 1 "reservoir"\n');
%  fprintf(fid,'2 2 "top"\n');
 fprintf(fid,'$EndPhysicalNames\n');
 fprintf(fid,'$Nodes\n');
 fprintf(fid,'%g\n',nTotNodes);
 
 for i = 1:nTotNodes
     fprintf(fid,'%g %g %g %g\n',i,Coords(i,1),Coords(i,2),Coords(i,3));
 end
 
 fprintf(fid,'$EndNodes\n');
 fprintf(fid,'$Elements\n');
 fprintf(fid,'%g\n',nTotEn);
 
if nElemNodes == 4
%   for j = 1:nTotSurfs
%      fprintf(fid,'%g %g %g %g %g %g %g %g\n',j,codeSup,surf_flag(j),surfTag,surf_PhyName,surfTopol(j,:));
%   end
  
  for m = 1:nTotElems
     fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',m,codeEl,elem_flag(m),cellTag,el_PhyName,elemTopol(m,:));
  end
end

if nElemNodes == 8
%   for j = 1:nTotSurfs
%      fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',j,codeSup,surf_flag(j),surfTag,surf_PhyName,surfTopol(j,:));
%   end
  
  for m = 1:nTotElems
     fprintf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g\n',m,codeEl,elem_flag,cellTag,el_PhyName,elemTopol(m,:));
  end
end
  
 fprintf(fid,'$EndElements\n');

 fclose(fid);