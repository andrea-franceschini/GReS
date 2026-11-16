function chronos(var)

   % Check if the variable is a char, if so change it to string
   if ischar(var)
      var = string(var);
   end

   % Check for the correct inputs
   if isstring(var) && ismember(var,["on","off","force"])
      use_Chronos = var;
   else
      error('Chronos expects a string which can be either "on", "off" or "force"');
   end

   % Assign value to a global variable
   setappdata(0,'use_Chronos', use_Chronos);
end