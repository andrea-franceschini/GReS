function varargout = pointToSurfaceProjection(P,n,coord)
      % Project nodes of triangle pair into auxiliary plane
      % get 2D direction of plane
      % Choose arbitrary vector not parallel to n
      if abs(n(1)) < 0.9
        temp = [1; 0; 0];
      else
        temp = [0; 1; 0];
      end
      % First direction in the plane
      d1 = cross(n, temp);
      d1 = d1 / norm(d1);
      % Second direction in the plane
      d2 = cross(n, d1);

      c = coord;
      cn = (c - P)*n;
      nN = size(c,1);
      projC = c - repmat(n',nN,1).*cn;
      % project in 2D
      x = (projC - P)*[d1 d2];

      varargout{1} = x;

      if nargout > 1
        varargout{2} = d1;
        varargout{3} = d2;
      end

    end