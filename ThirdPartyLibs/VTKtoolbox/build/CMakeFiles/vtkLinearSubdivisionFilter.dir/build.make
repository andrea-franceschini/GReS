# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build

# Include any dependencies generated for this target.
include CMakeFiles/vtkLinearSubdivisionFilter.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vtkLinearSubdivisionFilter.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vtkLinearSubdivisionFilter.dir/flags.make

CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o: CMakeFiles/vtkLinearSubdivisionFilter.dir/flags.make
CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o: ../vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o -c /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx

CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx > CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.i

CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx -o CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.s

CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o: CMakeFiles/vtkLinearSubdivisionFilter.dir/flags.make
CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o: /usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o   -c /usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c

CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c > CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.i

CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c -o CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.s

# Object files for target vtkLinearSubdivisionFilter
vtkLinearSubdivisionFilter_OBJECTS = \
"CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o" \
"CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o"

# External object files for target vtkLinearSubdivisionFilter
vtkLinearSubdivisionFilter_EXTERNAL_OBJECTS =

../MATLAB/vtkLinearSubdivisionFilter.mexa64: CMakeFiles/vtkLinearSubdivisionFilter.dir/vtkLinearSubdivisionFilter/vtkLinearSubdivisionFilter.cxx.o
../MATLAB/vtkLinearSubdivisionFilter.mexa64: CMakeFiles/vtkLinearSubdivisionFilter.dir/usr/local/MATLAB/R2023a/extern/version/c_mexapi_version.c.o
../MATLAB/vtkLinearSubdivisionFilter.mexa64: CMakeFiles/vtkLinearSubdivisionFilter.dir/build.make
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /usr/local/MATLAB/R2023a/extern/bin/glnxa64/libMatlabEngine.so
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /usr/local/MATLAB/R2023a/bin/glnxa64/libMatlabDataArray.so
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /usr/local/MATLAB/R2023a/bin/glnxa64/libmex.so
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /usr/local/MATLAB/R2023a/bin/glnxa64/libmx.so
../MATLAB/vtkLinearSubdivisionFilter.mexa64: libvtkMatlab.a
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersFlowPaths-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersParallel-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOGeometry-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOPLY-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersExtraction-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersModeling-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersHyperTree-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersTexture-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkParallelCore-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersHybrid-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersSources-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersGeneral-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersVerdict-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersGeometry-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkFiltersCore-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonComputationalGeometry-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOLegacy-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOXML-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOXMLParser-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkIOCore-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonExecutionModel-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonDataModel-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonTransforms-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonSystem-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonMisc-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonMath-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkCommonCore-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtktoken-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtkkissfft-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: /tmp/vtk/build/lib/libvtksys-9.3.so.9.3
../MATLAB/vtkLinearSubdivisionFilter.mexa64: CMakeFiles/vtkLinearSubdivisionFilter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library ../MATLAB/vtkLinearSubdivisionFilter.mexa64"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vtkLinearSubdivisionFilter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vtkLinearSubdivisionFilter.dir/build: ../MATLAB/vtkLinearSubdivisionFilter.mexa64

.PHONY : CMakeFiles/vtkLinearSubdivisionFilter.dir/build

CMakeFiles/vtkLinearSubdivisionFilter.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vtkLinearSubdivisionFilter.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vtkLinearSubdivisionFilter.dir/clean

CMakeFiles/vtkLinearSubdivisionFilter.dir/depend:
	cd /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles/vtkLinearSubdivisionFilter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vtkLinearSubdivisionFilter.dir/depend

