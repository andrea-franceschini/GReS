MATLAB_ROOT	= /usr/local/MATLAB/R2023a

MEX		= $(MATLAB_ROOT)/bin/mex

SOURCES		= mxImportGMSHmesh.cpp readGMSHmesh.cpp

HEADERS		= readGMSHmesh.hpp

MATLAB_FUNC = -DMEX_FUNCTION

all:
	@$(MEX) $(SOURCES) \
	CXXFLAGS="-I. -std=c++11 -O2 -fPIC $(MATLAB_FUNC)" \
	LDFLAGS="-L. -std=c++11 -O2"

clean:
	@rm -f mxImportGMSHmesh.mexa64
