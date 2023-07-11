function varargout = interpTable(xpt,ypt,dypt,x)
  varargout{1} = zeros(length(x),1);
  %
  x = -x;
  xNeg = x <= 0;
  xOutL = x > 0 & x < xpt(1);
  xIn = x >= xpt(1) & x <= xpt(end);
  xOutR = x > xpt(end);
  %
  varargout{1}(xNeg) = 1;
  %
  dy = (ypt(1) - 1)/xpt(1);
  varargout{1}(xOutL) = 1 + dy*x(xOutL);
  %
  varargout{1}(xIn) = interp1(xpt,ypt,x(xIn));
  [~,~,binID] = histcounts(x(xIn),xpt);
  %
  dy = (ypt(end) - ypt(end-1))/(xpt(end) - xpt(end-1));
  varargout{1}(xOutR) = ypt(end-1) + dy*(x(xOutR) - xpt(end-1));
  idOutR = find(xOutR);
  yNeg = varargout{1}(xOutR) < 0;
  varargout{1}(idOutR(yNeg)) = 0;
  if nargout == 2
    varargout{2} = zeros(length(x),1);
    %
    varargout{2}(xNeg) = 0;
    %
    varargout{2}(xOutL) = dy;
    %
    varargout{2}(xIn) = dypt(binID);
    %
    varargout{2}(idOutR) = dy;
    varargout{2}(idOutR(yNeg)) = 0;
  end
end

