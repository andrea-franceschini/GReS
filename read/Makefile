MATLAB_ROOT	= /usr/local/MATLAB/R2020b

MEX		= $(MATLAB_ROOT)/bin/mex

SOURCES		= mxImportGMSHmesh.cpp readGMSHmesh.cpp

HEADERS		= readGMSHmesh.hpp

all:
	@$(MEX) $(SOURCES) \
	CXXFLAGS="-I. -std=c++11 -O2 -fPIC" \
	LDFLAGS="-L. -std=c++11 -O2"

clean:
	@rm -f mxImportGMSHmesh.mexa64
