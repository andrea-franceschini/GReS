cmake_minimum_required(VERSION 3.11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
project(vtkToolbox)

# Move MEX binaries to vtkToolbox/MATLAB
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

## Set up Matlab dependency
find_package(Matlab REQUIRED)
if(NOT Matlab_FOUND)
    message(FATAL_ERROR "Matlab could not be found. Try setting Matlab_ROOT_DIR.")
endif()

## Set up VTK dependency
# Select only VTK components (libraries to be linked to) that are needed.
# (Removing the COMPONENTS argument from find_package links all of them.)
# See folder structure in https://vtk.org/doc/nightly/html/files.html
# to find out to which component a VTK class belongs.
set(VTK_COMPONENTS
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMath
    vtkCommonMisc
    vtkCommonSystem
    vtkFiltersCore
    vtkFiltersExtraction
    vtkFiltersFlowPaths
    vtkFiltersGeneral
    vtkFiltersGeometry
    vtkFiltersModeling
    vtkFiltersParallel
    vtkFiltersVerdict
    vtkIOCore
    vtkIOGeometry
    vtkIOLegacy
    vtkIOPLY
    vtkIOXML
)


find_package(VTK REQUIRED COMPONENTS ${VTK_COMPONENTS} NO_MODULE)
if(NOT VTK_FOUND)
    message(FATAL_ERROR "VTK could not be found. Try setting VTK_DIR.")
else()
    message(STATUS "Found VTK: " ${VTK_DIR})
endif()
if(VTK_VERSION LESS 8.9)
	include(${VTK_USE_FILE})
endif()

## Build library "vtkMatlab", which is used by all MEX functions in this toolbox
include_directories(libvtkMatlab ${Matlab_INCLUDE_DIRS})
add_library(vtkMatlab STATIC libvtkMatlab/common.cxx libvtkMatlab/vtkToStruct.cxx)
target_link_libraries(vtkMatlab ${VTK_LIBRARIES})

## Build MEX function vtkRead
matlab_add_mex(NAME vtkRead SRC vtkRead.cxx LINK_TO vtkMatlab ${VTK_LIBRARIES})