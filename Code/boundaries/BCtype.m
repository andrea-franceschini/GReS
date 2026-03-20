classdef BCtype

  % general class for handling a generic field in the problem

  enumeration
    dirichlet
    neumann
    source
  end

  methods (Static)

    % function out = isEssential(bctype)
    % 
    %   switch bctype
    %     case BCtype.dirichlet
    %       out = true;
    %     otherwise
    %       out = false;
    %   end
    % 
    %   if BCtype.isCustomBC(bctype)
    %     gresLog().warning(1, sprintf( ...
    %       "Custom boundary condition '%s' has been automatically considered to be NOT essential!\n" + ...
    %       "To override this behavior, specify the logical field 'essential' in the Boundary condition input.", ...
    %       bctype));
    %     out = false;
    %   end
    % 
    % end

    function out = isCustomBC(type)
      % detect if a Boundary condition is not a standard boundary condition
      try
        BCtype(type);
        out = false;
      catch
        out = true;
      end
    end


  end

end
