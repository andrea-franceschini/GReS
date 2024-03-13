#include "structToVtk.h"
#include "vtkToStruct.h"

#include <vtkCleanPolyData.h>

/* MATLAB entry function
 * nlhs/nrhs contain the number of left/right-hand-side arguments to this function
 * plhs/prhs are arrays of pointers to the arguments in MATLAB data format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ///// Parse inputs /////
    
    double tolerance = 0.0;
    bool toleranceIsAbsolute = false;
    
    const std::string syntax = "outStruct = vtkCleanPolyData(inStruct, (tolerance), (toleranceIsAbsolute))";
    if(nrhs < 1)
        mexErrMsgTxt(("Not enough input arguments. Syntax: " + syntax).c_str());
    if(nrhs > 1)
        tolerance = mxGetScalar(prhs[1]);
    if(nrhs > 2)
        toleranceIsAbsolute = mxGetScalar(prhs[2]);
    if(nrhs > 3)
        mexErrMsgTxt(("Too many input arguments. Syntax: " + syntax).c_str());
    if(nlhs > 1)
        mexErrMsgTxt(("Too many output arguments. Syntax: " + syntax).c_str());

    ///// Convert MATLAB struct into vtkPointSet /////
    
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkPointSet* pointSet = structToVtk(prhs[0], polyData, unstructuredGrid, false);
    
    ///// Apply vtkCleanPolyData /////
    
    if(pointSet->GetDataObjectType() != VTK_POLY_DATA)
        mexErrMsgTxt("vtkCleanPolyData requires poly data as input. Incompatible cell(s) found.");
    
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(pointSet);
    if(toleranceIsAbsolute)
    {
        cleanPolyData->ToleranceIsAbsoluteOn();
        cleanPolyData->SetAbsoluteTolerance(tolerance);
    }
    else
        cleanPolyData->SetTolerance(tolerance);
    cleanPolyData->Update();
    
    ///// Convert vtkPointSet into MATLAB struct /////
    
    plhs[0] = vtkToStruct(cleanPolyData->GetOutput(), false);
}
