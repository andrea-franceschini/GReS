MATLAB_ROOT	= /usr/local/MATLAB/R2020b

MEX		= $(MATLAB_ROOT)/bin/mex

SOURCES		= mxVTKWriter.cpp VTUWriter.cpp

HEADERS		= VTUWriter.hpp

all:
	@$(MEX) $(SOURCES) \
	CXXFLAGS="-I. -std=c++11 -O2 -fPIC" \
	LDFLAGS="-L. -std=c++11 -O2" \


clean:
	@rm -f mxVTKWriter.mexa64
