classdef Hexahedron < handle
  % Hexahedron element class

  properties (Access = public)
    % INPUT 
    dofmax = 3;
    % Number of faces of the element
    nFace = 6;
    % Number of nodes of the element
    nNode = [];
    % Node coordinates
    nodeCoords = [];
    % Total number of elements in the mesh
    nElem = [];
    % Element's topology
    elemTopol = [];
    % Cell centroid
    cellCentroid = [];
    % Surface's topology
    surfTopol = [];

    % OUTPUT 
    % Element side lenght
    Slenght = [];
    % Element volume
    Vol = [];
    % Element Area
    Area = [];
    % Shape functions matrix
    N = [];
    % Derivatives matrix
    B = [];
  end

  methods (Access = public)
    % Class constructor method
    function obj = Hexahedron(data)
       obj.setElementData(data);
    end
    
    function getSide(obj)
        % Allocation
        obj.Slenght = zeros(obj.nElem,3);
        for el = 1 : obj.nElem %element row in element topology
              % Inizialyzing variables
              a = 0;
              b = 0;
              c = 0;
              j = obj.elemTopol(el,1); % first node id in element topology
              for i = 2 :obj.nNode
                 m = obj.elemTopol(el,i); %cycle for all element node ids in element topology
                 if (obj.nodeCoords(j,1)) ~= (obj.nodeCoords(m,1)) && a==0
                 a = abs(obj.nodeCoords(m,1)-obj.nodeCoords(j,1));
                 elseif (obj.nodeCoords(j,2)) ~= (obj.nodeCoords(m,2)) && b==0
                 b = abs(obj.nodeCoords(m,2)-obj.nodeCoords(j,2));
                 elseif (obj.nodeCoords(j,3)) ~= (obj.nodeCoords(m,3)) && c==0
                 c = abs(obj.nodeCoords(m,3)-obj.nodeCoords(j,3));
                 end
               end
              if a~=0 && b~=0 && c~=0
              obj.Slenght(el,:) = [a b c];
              else
                  error('Error computing element side length')
              end
        end
        
    end
    
    
    % Function to re order the element topology based on the natural coordinate 
    % system
    function obj = reOrderTopol(obj)
      old_topol = obj.elemTopol;
      new_topol = zeros(obj.nElem,obj.nNode);
      for el = 1:obj.nElem
          for i = 1:obj.nNode
              node = old_topol(el,i);
              x_coords(1,i) = obj.nodeCoords(node,1);
              y_coords(1,i) = obj.nodeCoords(node,2);
              z_coords(1,i) = obj.nodeCoords(node,3);
          end
          x_min = min(x_coords);
          y_min = min(y_coords);
          z_min = min(z_coords);
          
          for i = 1:obj.nNode
              node = old_topol(el,i);
              c1 = obj.nodeCoords(node,1);
              c2 = obj.nodeCoords(node,2);
              c3 = obj.nodeCoords(node,3);
              if c1 == x_min && c2 == y_min && c3 == z_min
                 new_topol(el,1) = node;
              elseif c1 ~= x_min && c2 == y_min && c3 == z_min
                 new_topol(el,2) = node;
              elseif c1 ~= x_min && c2 ~= y_min && c3 == z_min
                 new_topol(el,3) = node;
              elseif c1 == x_min && c2 ~= y_min && c3 == z_min
                 new_topol(el,4) = node;
              elseif c1 == x_min && c2 == y_min && c3 ~= z_min
                 new_topol(el,5) = node;
              elseif c1 ~= x_min && c2 == y_min && c3 ~= z_min
                 new_topol(el,6) = node;
              elseif c1 ~= x_min && c2 ~= y_min && c3 ~= z_min
                 new_topol(el,7) = node;
              elseif c1 == x_min && c2 ~= y_min && c3 ~= z_min
                 new_topol(el,8) = node;
              end
          end
          obj.elemTopol(el,:) = new_topol(el,:);
      end
    end
    
    

    %Function for shape function coefficients calculation
    function getCoefficient(obj)
        % Allocation
        obj.N = zeros(obj.dofmax*obj.nElem,obj.dofmax*obj.nNode);
        
        % Linear hexahedron
        if obj.nNode == 8
            
          for el = 1:obj.nElem
            for j = 1:obj.nNode
              nod = obj.elemTopol(el,j);
              
              s = (obj.nodeCoords(nod,1)-obj.cellCentroid(el,1))/(0.5*obj.Slenght(el,1));
              t = (obj.nodeCoords(nod,2)-obj.cellCentroid(el,2))/(0.5*obj.Slenght(el,2));
              p = (obj.nodeCoords(nod,3)-obj.cellCentroid(el,3))/(0.5*obj.Slenght(el,3));
              
              N1 = 1/8*(1-s)*(1-t)*(1-p);
              N2 = 1/8*(1+s)*(1-t)*(1-p);
              N3 = 1/8*(1+s)*(1+t)*(1-p);
              N4 = 1/8*(1-s)*(1+t)*(1-p);
              N5 = 1/8*(1-s)*(1-t)*(1+p);
              N6 = 1/8*(1+s)*(1-t)*(1+p);
              N7 = 1/8*(1+s)*(1+t)*(1+p);
              N8 = 1/8*(1-s)*(1+t)*(1+p);
              
              totN = N1+N2+N3+N4+N5+N6+N7+N8;
              if totN ~= 1
                  error('Shape functions must have local support')
              end
              
              obj.N(3*el-2:3*el,3*j-2:3*j) = [totN 0 0;
                                              0 totN 0;
                                              0 0 totN];
            end
          end


        % Quadratic hexahedron
        elseif obj.nNode == 20
             error('Element not available');
        end
    end
    
    
    %Function for volume calculation
    function getVolume(obj)
        %Allocation
        obj.Vol = zeros(obj.nElem,1);
        
        for el = 1:obj.nElem
            a = obj.Slenght(el,1);
            b = obj.Slenght(el,2);
            c = obj.Slenght(el,3);
            
            obj.Vol(el) = a*b*c;
         end
    end
    
    
       %Function for derivatives matrices calculation
     function getDerivatives(obj)
       %Allocation
       obj.B = zeros(2*obj.dofmax*obj.nElem,obj.nNode*obj.dofmax);
         
       % Linear hexahedron
       if obj.nNode == 8
            
          for el = 1:obj.nElem
            for j = 1:obj.nNode
              nod = obj.elemTopol(el,j);
              
              s = (obj.nodeCoords(nod,1)-obj.cellCentroid(el,1))/(0.5*obj.Slenght(el,1));
              t = (obj.nodeCoords(nod,2)-obj.cellCentroid(el,2))/(0.5*obj.Slenght(el,2));
              p = (obj.nodeCoords(nod,3)-obj.cellCentroid(el,3))/(0.5*obj.Slenght(el,3));
                         
              if j == 1 %first node of the element
              b = -1/(0.5*obj.Slenght(el,1))*(1-t)*(1-p);
              c = -1/(0.5*obj.Slenght(el,2))*(1-s)*(1-p);
              d = -1/(0.5*obj.Slenght(el,3))*(1-s)*(1-t); 
              
              elseif j == 2 %second node of the element
              b = 1/(0.5*obj.Slenght(el,1))*(1-t)*(1-p);
              c = -1/(0.5*obj.Slenght(el,2))*(1+s)*(1-p);
              d = -1/(0.5*obj.Slenght(el,3))*(1+s)*(1-t); 
              
              elseif j == 3 %third node of the element
              b = 1/(0.5*obj.Slenght(el,1))*(1+t)*(1-p);
              c = 1/(0.5*obj.Slenght(el,2))*(1+s)*(1-p);
              d = -1/(0.5*obj.Slenght(el,3))*(1+s)*(1+t); 
              
              elseif j == 4 %fourth node of the element
              b = -1/(0.5*obj.Slenght(el,1))*(1+t)*(1-p);
              c = 1/(0.5*obj.Slenght(el,2))*(1-s)*(1-p);
              d = -1/(0.5*obj.Slenght(el,3))*(1-s)*(1+t); 
              
              elseif j == 5 %fifth node of the element
              b = -1/(0.5*obj.Slenght(el,1))*(1-t)*(1+p);
              c = -1/(0.5*obj.Slenght(el,2))*(1-s)*(1+p);
              d = 1/(0.5*obj.Slenght(el,3))*(1-s)*(1-t); 
              
              elseif j == 6 %sixth node of the element
              b = 1/(0.5*obj.Slenght(el,1))*(1-t)*(1+p);
              c = -1/(0.5*obj.Slenght(el,2))*(1+s)*(1+p);
              d = 1/(0.5*obj.Slenght(el,3))*(1+s)*(1-t); 
              
              elseif j == 7 %seventh node of the element
              b = 1/(0.5*obj.Slenght(el,1))*(1+t)*(1+p);
              c = 1/(0.5*obj.Slenght(el,2))*(1+s)*(1+p);
              d = 1/(0.5*obj.Slenght(el,3))*(1+s)*(1+t); 
              
              elseif j == 8 %eighth node of the element
              b = -1/(0.5*obj.Slenght(el,1))*(1+t)*(1+p);
              c = 1/(0.5*obj.Slenght(el,2))*(1-s)*(1+p);
              d = 1/(0.5*obj.Slenght(el,3))*(1-s)*(1+t); 
              end
              
              B = [b    0    0;
                   0    c    0;
                   0    0    d;
                   c    b    0;
                   d    0    b;
                   0    d    c];
               
              obj.B (6*el-5:6*el,3*j-2:3*j)  = 1/8 * B; 
              
            end
          end

         
       % Quadratic hexahedron
       elseif obj.nNode == 20
             error('Element not available');
       end
      end
       
  
  % Function for the surfaces area calculation
    function getArea(obj)
    [nrow,ncol] = size(obj.surfTopol);
     obj.Area = zeros(nrow,1);
         for s = 1:nrow
             for i = 1:ncol
                 node = obj.surfTopol(s,i);
                 coords(i,:) = obj.nodeCoords(node,:);
             end 

             L = zeros (1,obj.dofmax);
             for k = 1:obj.dofmax
                 for i = 1:ncol
                     if i == ncol
                        break
                     elseif coords(i,k) ~= coords(i+1,k)
                        L(k) = abs(coords(i,k) - coords(i+1,k));
                        if L(k) ~= 0
                           break
                       end
                     end
                 end
             end
             L(L==0) = [];
             obj.Area(s) = L(1)*L(2);
         end
    end
    
  end


  methods (Access = private)
      % Function that set the element parameters coming from "data"
    function setElementData(obj,data)
        nNode = data{1};
        obj.nNode = nNode(1,1);
        obj.nodeCoords = data{2};
        obj.nElem = data{3};
        obj.elemTopol = data{4};
        obj.surfTopol = data{5};
        obj.cellCentroid = data{6};
    end
                
                
  end

end
