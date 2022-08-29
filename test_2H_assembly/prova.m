close all;
clear;
clc;


nodeCoords = [3 2;
             5 2;
             5 4;
             3 4];
s = -0.5773;
t = -0.5773;

nNode = 4;
for n = 1:nNode
                   if n == 1 %first node of the element
                   dNs = -1/4*(1-t);
                   dNt = -1/4*(1-s);
                 
              
                   elseif n == 2 %second node of the element
                   dNs = 1/4*(1-t);
                   dNt = -1/4*(1+s);
                   
              
                   elseif n == 3 %third node of the element
                   dNs = 1/4*(1+t);
                   dNt = 1/4*(1+s);
                    
              
                   elseif n == 4 %fourth node of the element
                   dNs = -1/4*(1+t);
                   dNt = 1/4*(1-s);
                    
              
                   end
          dN = [dNs; dNt];
          derMatrix (:,n) = dN;
          coordMatrix(n,:) = nodeCoords(n,:);
end
            
J = derMatrix*coordMatrix;
detJ = det(J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:nNode 
              
              invJ = inv(J);
              
              if n == 1 %first node of the element
              b = invJ(1,:)*derMatrix(:,n);
              c = invJ(2,:)*derMatrix(:,n);
              
              elseif n == 2 %second node of the element
              b = invJ(1,:)*derMatrix(:,n);
              c = invJ(2,:)*derMatrix(:,n); 
              
              elseif n == 3 %third node of the element
              b = invJ(1,:)*derMatrix(:,n);
              c = invJ(2,:)*derMatrix(:,n); 
              
              elseif n == 4 %fourth node of the element
              b = invJ(1,:)*derMatrix(:,n);
              c = invJ(2,:)*derMatrix(:,n); 
               
              end             
              
              B(:,2*n-1:2*n) = [b    0;
                                0    c;
                                c    b];
end
               
       B = 1/detJ * B; 
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa = 1/4*(nodeCoords(1,2)*(s-1)+nodeCoords(2,2)*(-1-s)+nodeCoords(3,2)*(1+s)+nodeCoords(4,2)*(1-s));
bb = 1/4*(nodeCoords(1,2)*(t-1)+nodeCoords(2,2)*(1-t)+nodeCoords(3,2)*(1+t)+nodeCoords(4,2)*(-1-t));
cc = 1/4*(nodeCoords(1,1)*(t-1)+nodeCoords(2,1)*(1-t)+nodeCoords(3,1)*(1+t)+nodeCoords(4,1)*(-1-t));
dd = 1/4*(nodeCoords(1,1)*(s-1)+nodeCoords(2,1)*(-1-s)+nodeCoords(3,1)*(1+s)+nodeCoords(4,1)*(1-s));


N1s = -0.394325;
N1t = -0.394325;

N2s = 1/4*(1-t);
N2t = -1/4*(1+s);
N3s = 1/4*(1+t);
N3t = 1/4*(1+s);
N4s = -1/4*(1+t);
N4t = 1/4*(1-s);

B1 = [aa*N1s-bb*N1t 0;
      0 cc*N1t-dd*N1s;
      cc*N1t-dd*N1s aa*N1s-bb*N1t];
B2 = [aa*N2s-bb*N2t 0;
      0 cc*N2t-dd*N2s;
      cc*N2t-dd*N2s aa*N2s-bb*N2t];
B3 = [aa*N3s-bb*N3t 0;
      0 cc*N3t-dd*N3s;
      cc*N3t-dd*N3s aa*N3s-bb*N3t];
B4 = [aa*N4s-bb*N4t 0;
      0 cc*N4t-dd*N4s;
      cc*N4t-dd*N4s aa*N4s-bb*N4t];
       