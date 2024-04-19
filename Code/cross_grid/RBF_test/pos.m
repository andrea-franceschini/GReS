function val = pos(val)
    %input is an array; convert to 0 all negative values
    val(val<0) = 0; 
end

