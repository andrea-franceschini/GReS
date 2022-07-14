classdef Tetrahedron < handle
  % Tetrahedron element class

  properties (Access = public)
    % INPUT 
    dofmax = 3;
    % Number of faces of the element
    nFace = 4;
    % Number of nodes of the element
    nNode = [];
    % Node coordinates
    nodeCoords = [];
    % Total number of elements in the mesh
    nElem = [];
    % Element's topology
    elemTopol = [];
    % Surface's topology
    surfTopol = [];

    % OUTPUT 
    % Element volume
    Vol = [];
    % Element Area
    Area = [];
    % Shape function coefficients
    Coeff = [];
    % Derivates matrices
    B = [];
  end

  methods (Access = public)
    % Class constructor method
    function obj = Tetrahedron(data)
       obj.setElementData(data);
    end

    
    %Function for shape function coefficients calculation
    function getCoefficient(obj)
        % Allocation
        obj.Coeff = zeros(3*obj.nElem,obj.nNode);
        
        % Linear tetrahedron
        if obj.nNode == 4
            
          for el = 1:obj.nElem
            i = obj.elemTopol(el,1);
            j = obj.elemTopol(el,2);
            m = obj.elemTopol(el,3);
            p = obj.elemTopol(el,4);
            
            
            b(1) = -det([1 obj.nodeCoords(j,2) obj.nodeCoords(j,3);
                         1 obj.nodeCoords(m,2) obj.nodeCoords(m,3);
                         1 obj.nodeCoords(p,2) obj.nodeCoords(p,3);]);
            c(1) = det([1 obj.nodeCoords(j,1) obj.nodeCoords(j,3);
                        1 obj.nodeCoords(m,1) obj.nodeCoords(m,3);
                        1 obj.nodeCoords(p,1) obj.nodeCoords(p,3);]);
            d(1) = -det([1 obj.nodeCoords(j,1) obj.nodeCoords(j,2);
                         1 obj.nodeCoords(m,1) obj.nodeCoords(m,2);
                         1 obj.nodeCoords(p,1) obj.nodeCoords(p,2);]);
                     
                     
            b(2) = det([1 obj.nodeCoords(i,2) obj.nodeCoords(i,3);
                        1 obj.nodeCoords(m,2) obj.nodeCoords(m,3);
                        1 obj.nodeCoords(p,2) obj.nodeCoords(p,3);]);
            c(2) = -det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,3);
                         1 obj.nodeCoords(m,1) obj.nodeCoords(m,3);
                         1 obj.nodeCoords(p,1) obj.nodeCoords(p,3);]);
            d(2) = det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,2);
                        1 obj.nodeCoords(m,1) obj.nodeCoords(m,2);
                        1 obj.nodeCoords(p,1) obj.nodeCoords(p,2);]);
                                         
                     
            b(3) = -det([1 obj.nodeCoords(i,2) obj.nodeCoords(i,3);
                         1 obj.nodeCoords(j,2) obj.nodeCoords(j,3);
                         1 obj.nodeCoords(p,2) obj.nodeCoords(p,3);]);
            c(3) = det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,3);
                        1 obj.nodeCoords(j,1) obj.nodeCoords(j,3);
                        1 obj.nodeCoords(p,1) obj.nodeCoords(p,3);]);
            d(3) = -det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,2);
                         1 obj.nodeCoords(j,1) obj.nodeCoords(j,2);
                         1 obj.nodeCoords(p,1) obj.nodeCoords(p,2);]);
                     
                     
            b(4) = det([1 obj.nodeCoords(i,2) obj.nodeCoords(i,3);
                        1 obj.nodeCoords(j,2) obj.nodeCoords(j,3);
                        1 obj.nodeCoords(m,2) obj.nodeCoords(m,3);]);
            c(4) = -det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,3);
                         1 obj.nodeCoords(j,1) obj.nodeCoords(j,3);
                         1 obj.nodeCoords(m,1) obj.nodeCoords(m,3);]);
            d(4) = det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,2);
                        1 obj.nodeCoords(j,1) obj.nodeCoords(j,2);
                        1 obj.nodeCoords(m,1) obj.nodeCoords(m,2);]);
                                         
            obj.Coeff(3*el-2,:) = b;
            obj.Coeff(3*el-1,:) = c;
            obj.Coeff(3*el,:) = d;
          end

        % Quadratic tetrahedron
        elseif obj.nNode == 10
             error('Element not available');
        end
    end
    
    
    %Function for volume calculation
    function getVolume(obj)
        
        for el = 1:obj.nElem
            i = obj.elemTopol(el,1);
            j = obj.elemTopol(el,2);
            m = obj.elemTopol(el,3);
            p = obj.elemTopol(el,4);
            
            obj.Vol = (det([1 obj.nodeCoords(i,1) obj.nodeCoords(i,2) obj.nodeCoords(i,3);
                        1 obj.nodeCoords(j,1) obj.nodeCoords(j,2) obj.nodeCoords(j,3);
                        1 obj.nodeCoords(m,1) obj.nodeCoords(m,2) obj.nodeCoords(m,3);
                        1 obj.nodeCoords(p,1) obj.nodeCoords(p,2) obj.nodeCoords(p,3)]))/6;
         end
    end
    
    
        %Function for derivates matrices calculation
     function getDerivates(obj)
        %Allocation
        obj.B = zeros(2*obj.dofmax*obj.nElem,obj.nNode*obj.dofmax);
         
       % Linear tetrahedron
       if obj.nNode == 4
            
         for el = 1: obj.nElem
          b = obj.Coeff(3*el-2,:);
          c = obj.Coeff(3*el-1,:);
          d = obj.Coeff(3*el,:);

             B1 = [b(1)   0     0;
                   0    c(1)    0;
                   0    0    d(1);
                   c(1) b(1)    0;
                   0    d(1) c(1);
                   d(1) 0    b(1)];
                          
             B2 = [b(2)   0     0;
                   0    c(2)    0;
                   0    0    d(2);
                   c(2) b(2)    0;
                   0    d(2) c(2);
                   d(2) 0    b(2)];
                          
             B3 = [b(3)   0     0;
                   0    c(3)    0;
                   0    0    d(3);
                   c(3) b(3)    0;
                   0    d(3) c(3);
                   d(3) 0    b(3)];
               
             B4 = [b(4)   0     0;
                   0    c(4)    0;
                   0    0    d(4);
                   c(4) b(4)    0;
                   0    d(4) c(4);
                   d(4) 0    b(4)];
                          
         obj.B (6*el-5:6*el,:)  = (1/(6*obj.Vol))* [B1 B2 B3 B4];
         end
         
       % Quadratic tetrahedron
       elseif obj.nNode == 10
             error('Element not available');
       end
      end
       
  
  % Function for the surfaces area calculation
    function getArea(obj)
    [nrow,ncol] = size(obj.surfTopol);
     obj.Area = zeros(nrow,1);
         for f = 1:nrow
             for ind = 1:ncol
                 node = obj.surfTopol(f,ind);
                 coord(ind,:) = obj.nodeCoords(node,:);
             end
             n1 = coord(1,:);
             n2 = coord(2,:);
             n3 = coord(3,:);

             v1 = [n2(1)-n1(1),n2(2)-n1(2),n2(3)-n1(3)];
             v2 = [n3(1)-n1(1),n3(2)-n1(2),n3(3)-n1(3)];
             
             prodvett = cross(v1,v2);
             norma = norm(prodvett);
             obj.Area(f) = 0.5*norma;
            
         end
     end
        
     end


  methods (Access = private)
      % Function that set the element parameters coming from "data"
    function setElementData(obj,data)
        obj.nNode = data{1};
        obj.nodeCoords = data{2};
        obj.nElem = data{3};
        obj.elemTopol = data{4};
        obj.surfTopol = data{5};
    end
  end

end
