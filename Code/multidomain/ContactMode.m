classdef ContactMode
  
  enumeration
    stick       % 1
    newSlip     % 2
    slip        % 3
    open        % 4
  end

  methods (Static)
    function intModes = integer(modes)
      % convert the enumeration into integers 
      intModes = zeros(numel(modes),1);
      intModes(modes == ContactMode.stick) = 1;
      intModes(modes == ContactMode.newSlip) = 2;
      intModes(modes == ContactMode.slip) = 3;
      intModes(modes == ContactMode.open) = 4;
    end
  end


end