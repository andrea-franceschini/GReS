function [newTopol] = reOrderTopol(oldTopol,coord,nElem)
      nNode = 8;
      newTopol = zeros(nElem,nNode);
      for el = 1:nElem
          x_coords = zeros(nNode,1);
          y_coords = zeros(nNode,1);
          z_coords = zeros(nNode,1);
          for i = 1:nNode
             node = oldTopol(el,i);
             x_coords(i) = coord(node,1);
             y_coords(i) = coord(node,2);
             z_coords(i) = coord(node,3);
          end
          x_min = min(x_coords);
          y_min = min(y_coords);
          z_min = min(z_coords);
          
          for i = 1:nNode
              node = oldTopol(el,i);
              c1 = coord(node,1);
              c2 = coord(node,2);
              c3 = coord(node,3);
              if c1 == x_min && c2 == y_min && c3 == z_min
                 newTopol(el,1) = node;
              elseif c1 ~= x_min && c2 == y_min && c3 == z_min
                 newTopol(el,2) = node;
              elseif c1 ~= x_min && c2 ~= y_min && c3 == z_min
                 newTopol(el,3) = node;
              elseif c1 == x_min && c2 ~= y_min && c3 == z_min
                 newTopol(el,4) = node;
              elseif c1 == x_min && c2 == y_min && c3 ~= z_min
                 newTopol(el,5) = node;
              elseif c1 ~= x_min && c2 == y_min && c3 ~= z_min
                 newTopol(el,6) = node;
              elseif c1 ~= x_min && c2 ~= y_min && c3 ~= z_min
                 newTopol(el,7) = node;
              elseif c1 == x_min && c2 ~= y_min && c3 ~= z_min
                 newTopol(el,8) = node;
              end
          end
%           newTopol(el,:) = newTopol(el,:);
      end
end