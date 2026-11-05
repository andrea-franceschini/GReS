cd mex
mex -R2018a -O CXXFLAGS="$CXXFLAGS -std=c++14" mxVTKWriter.cpp VTUWriter.cpp