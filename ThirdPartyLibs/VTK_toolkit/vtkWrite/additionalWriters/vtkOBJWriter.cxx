////////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2013, OpenCV Foundation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
// Author: Anatoly Baksheev, Itseez Inc.
//

#include "vtkOBJWriter.h"

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkInformation.h>
#include <vtkErrorCode.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkVersionMacros.h>

vtkStandardNewMacro(vtkOBJWriter);

vtkOBJWriter::vtkOBJWriter()
{
    std::ofstream fout; // only used to extract the default precision
    this->DecimalPrecision = fout.precision();
    this->FileName = NULL;
}

vtkOBJWriter::~vtkOBJWriter(){}

void vtkOBJWriter::WriteData()
{
    vtkPolyData *input = this->GetInput();
    if(!input)
        return;
    
    if(!this->FileName)
    {
        vtkErrorMacro(<< "No FileName specified! Can't write!");
        this->SetErrorCode(vtkErrorCode::NoFileNameError);
        return;
    }
    
    vtkDebugMacro(<<"Opening vtk file for writing...");
    ostream *outfilep = new std::ofstream(this->FileName, ios::out);
    if(outfilep->fail())
    {
        vtkErrorMacro(<< "Unable to open file: "<< this->FileName);
        this->SetErrorCode(vtkErrorCode::CannotOpenFileError);
        delete outfilep;
        return;
    }
    
    std::ostream& outfile = *outfilep;
    double p[3];
    
    // write points
    for(int i = 0; i < input->GetNumberOfPoints(); i++)
    {
        input->GetPoint(i, p);
        outfile << std::setprecision(this->DecimalPrecision) << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
    
    const int idStart = 1;
    
    // write point data
    vtkSmartPointer<vtkDataArray> normals = input->GetPointData()->GetNormals();
    if(normals)
    {
        for(int i = 0; i < normals->GetNumberOfTuples(); i++)
        {
            normals->GetTuple(i, p);
            outfile << std::setprecision(this->DecimalPrecision) << "vn " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }
    
    vtkSmartPointer<vtkDataArray> tcoords = input->GetPointData()->GetTCoords();
    if(tcoords)
    {
        for(int i = 0; i < tcoords->GetNumberOfTuples(); i++)
        {
            tcoords->GetTuple(i, p);
            outfile << std::setprecision(this->DecimalPrecision) << "vt " << p[0] << " " << p[1] << std::endl;
        }
    }
    
    // write verts, if any
    if(input->GetNumberOfVerts() > 0)
    {
        vtkIdType npts = 0;
        #if VTK_MAJOR_VERSION > 8
            const vtkIdType *index = 0;
        #else
               vtkIdType *index = 0;
        #endif
        vtkCellArray *cells = input->GetVerts();
        for(cells->InitTraversal(); cells->GetNextCell(npts, index); )
        {
            outfile << "p ";
            for(int i = 0; i < npts; i++)
                outfile << index[i] + idStart << " ";
            outfile << std::endl;
        }
    }
    
    // write lines, if any
    if(input->GetNumberOfLines() > 0)
    {
        vtkIdType npts = 0;
        #if VTK_MAJOR_VERSION > 8
            const vtkIdType *index = 0;
        #else
               vtkIdType *index = 0;
        #endif
        vtkCellArray *cells = input->GetLines();
        for(cells->InitTraversal(); cells->GetNextCell(npts, index); )
        {
            outfile << "l ";
            if(tcoords)
            {
                for(int i = 0; i < npts; i++)
                    outfile << index[i] + idStart << "/" << index[i] + idStart << " ";
            }
            else
                for(int i = 0; i < npts; i++)
                    outfile << index[i] + idStart << " ";
            
            outfile << std::endl;
        }
    }
    
    // write polys, if any
    if(input->GetNumberOfPolys() > 0)
    {
        vtkIdType npts = 0;
        #if VTK_MAJOR_VERSION > 8
            const vtkIdType *index = 0;
        #else
               vtkIdType *index = 0;
        #endif
        vtkCellArray *cells = input->GetPolys();
        for(cells->InitTraversal(); cells->GetNextCell(npts, index); )
        {
            outfile << "f ";
            for(int i = 0; i < npts; i++)
            {
                if(normals)
                {
                    if(tcoords)
                        outfile << index[i] + idStart << "/"  << index[i] + idStart << "/" << index[i] + idStart << " ";
                    else
                        outfile << index[i] + idStart << "//" << index[i] + idStart << " ";
                }
                else
                {
                    if(tcoords)
                        outfile << index[i] + idStart << " " << index[i] + idStart << " ";
                    else
                        outfile << index[i] + idStart << " ";
                }
            }
            outfile << std::endl;
        }
    }
    
    // write tstrips, if any
    if(input->GetNumberOfStrips() > 0)
    {
        vtkIdType npts = 0;
        #if VTK_MAJOR_VERSION > 8
            const vtkIdType *index = 0;
        #else
               vtkIdType *index = 0;
        #endif
        vtkCellArray *cells = input->GetStrips();
        for(cells->InitTraversal(); cells->GetNextCell(npts, index); )
        {
            for(int i = 2, i1, i2; i < npts; ++i)
            {
                if(i % 2)
                {
                    i1 = i - 1;
                    i2 = i - 2;
                }
                else
                {
                    i1 = i - 1;
                    i2 = i - 2;
                }
                
                if(normals)
                {
                    if(tcoords)
                    {
                        outfile << "f " << index[i1] + idStart << "/" << index[i1] + idStart << "/" << index[i1] + idStart << " "
                        << index[i2]+ idStart << "/" << index[i2] + idStart << "/" << index[i2] + idStart << " "
                        << index[i] + idStart << "/" << index[i]  + idStart << "/" << index[i]  + idStart << std::endl;
                    }
                    else
                    {
                        outfile << "f " << index[i1] + idStart << "//" << index[i1] + idStart << " " << index[i2] + idStart
                        << "//" << index[i2] + idStart << " "  << index[i]  + idStart << "//" << index[i] + idStart << std::endl;
                    }
                }
                else
                {
                    if(tcoords)
                    {
                        outfile << "f " << index[i1] + idStart << "/" << index[i1] + idStart << " " << index[i2] + idStart
                        << "/" << index[i2] + idStart << " "  << index[i]  + idStart << "/" << index[i]  + idStart << std::endl;
                    }
                    else
                        outfile << "f " << index[i1] + idStart << " " << index[i2] + idStart << " " << index[i] + idStart << std::endl;
                }
            }
        }
    }
    
    vtkDebugMacro(<<"Closing vtk file\n");
    delete outfilep;
    
    // Delete the file if an error occurred
    if(this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
        vtkErrorMacro("Ran out of disk space; deleting file: " << this->FileName);
        std::remove(this->FileName);
    }
}

void vtkOBJWriter::PrintSelf(ostream& os, vtkIndent indent)
{
    Superclass::PrintSelf(os, indent);
    os << indent << "DecimalPrecision: " << DecimalPrecision << "\n";
}

int vtkOBJWriter::FillInputPortInformation(int, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}

vtkPolyData* vtkOBJWriter::GetInput()
{
    return vtkPolyData::SafeDownCast(this->Superclass::GetInput());
}

vtkPolyData* vtkOBJWriter::GetInput(int port)
{
    return vtkPolyData::SafeDownCast(this->Superclass::GetInput(port));
}
