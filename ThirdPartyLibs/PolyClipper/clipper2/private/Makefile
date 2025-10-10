# Change the next line to the MATLAB application
MATLAB=/Applications/MATLAB_R2016b.app

# Change the next line to the clipper directory (e.g., "clipper", "clipper_ver6", etc.)
CLIPPER_DIR=clipper



MEX=$(MATLAB)/bin/mex
FLAGS=-D__int64=__int64_t
# OLDFLAGS=-Dchar16_t=uint16_t -D__int64=__int64_t

clipper.mexmaci64: clipper.cpp clipper.hpp mexclipper.cpp
	$(MEX) $(FLAGS) clipper.cpp mexclipper.cpp

clipper.cpp: $(CLIPPER_DIR)/cpp/clipper.cpp
	cp $(CLIPPER_DIR)/cpp/clipper.cpp .

clipper.hpp: $(CLIPPER_DIR)/cpp/clipper.hpp
	cp $(CLIPPER_DIR)/cpp/clipper.hpp .
