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
include CMakeFiles/vtkMatlab.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vtkMatlab.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vtkMatlab.dir/flags.make

CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o: CMakeFiles/vtkMatlab.dir/flags.make
CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o: ../libvtkMatlab/common.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o -c /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/common.cxx

CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/common.cxx > CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.i

CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/common.cxx -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.s

CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o: CMakeFiles/vtkMatlab.dir/flags.make
CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o: ../libvtkMatlab/vtkToStruct.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o -c /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/vtkToStruct.cxx

CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/vtkToStruct.cxx > CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.i

CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/vtkToStruct.cxx -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.s

CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o: CMakeFiles/vtkMatlab.dir/flags.make
CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o: ../libvtkMatlab/structToVtk.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o -c /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/structToVtk.cxx

CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/structToVtk.cxx > CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.i

CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/libvtkMatlab/structToVtk.cxx -o CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.s

# Object files for target vtkMatlab
vtkMatlab_OBJECTS = \
"CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o" \
"CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o" \
"CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o"

# External object files for target vtkMatlab
vtkMatlab_EXTERNAL_OBJECTS =

libvtkMatlab.a: CMakeFiles/vtkMatlab.dir/libvtkMatlab/common.cxx.o
libvtkMatlab.a: CMakeFiles/vtkMatlab.dir/libvtkMatlab/vtkToStruct.cxx.o
libvtkMatlab.a: CMakeFiles/vtkMatlab.dir/libvtkMatlab/structToVtk.cxx.o
libvtkMatlab.a: CMakeFiles/vtkMatlab.dir/build.make
libvtkMatlab.a: CMakeFiles/vtkMatlab.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libvtkMatlab.a"
	$(CMAKE_COMMAND) -P CMakeFiles/vtkMatlab.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vtkMatlab.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vtkMatlab.dir/build: libvtkMatlab.a

.PHONY : CMakeFiles/vtkMatlab.dir/build

CMakeFiles/vtkMatlab.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vtkMatlab.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vtkMatlab.dir/clean

CMakeFiles/vtkMatlab.dir/depend:
	cd /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build /scratches/russel3/moretto/GReS/ThirdPartyLibs/VTKtoolbox/build/CMakeFiles/vtkMatlab.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vtkMatlab.dir/depend

