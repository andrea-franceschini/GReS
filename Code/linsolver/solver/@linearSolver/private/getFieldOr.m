function val = getFieldOr(s, f, defaultVal)
% GETFIELDOR Safely retrieves struct field or returns default value if absent or empty

   if isfield(s, f) && ~isempty(s.(f))
      val = s.(f);
   else
      val = defaultVal;
   end
end
