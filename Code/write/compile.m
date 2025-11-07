cd mex
fName = "mxVTKWriter";
if exist("mxVTKWriter","file") == 3
  mex -O mxVTKWriter.cpp VTUWriter.cpp
else
  f = fName+"mex";
  disp("Mex file %s already exists. Compilation is skipped")
end
