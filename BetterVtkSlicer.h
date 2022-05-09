//
// Created by peter on 09.05.2022.
//

#ifndef VIS2021_BETTERVTKSLICER_H
#define VIS2021_BETTERVTKSLICER_H

#include "vtkImageAlgorithm.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkInformationVector.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "utils.h"

class BetterVtkSlicer : public vtkImageAlgorithm {
public:
    static BetterVtkSlicer *New();

vtkTypeMacro(BetterVtkSlicer, vtkImageAlgorithm);

    BetterVtkSlicer() {}

    void SetHeight(int h) {
        height = h;
    }

    void SetZSlice() {
        zSlice = true;
        xSlice = false;
        ySlice = false;
    }

    void SetXSlice() {
        zSlice = false;
        xSlice = true;
        ySlice = false;
    }

    void SetYSlice() {
        ySlice = true;
        xSlice = false;
        zSlice = false;
    }


protected:
    int RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector **inputVector,
                    vtkInformationVector *outputVector) {
        vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
        vtkInformation *outInfo = outputVector->GetInformationObject(0);
        vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
        vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
        vtkPointData *inPointData = input->GetPointData();
        vtkPointData *outPointData = output->GetPointData();

        int dimensions[3];
        int newdimensions[3];
        input->GetDimensions(dimensions);
        input->GetDimensions(newdimensions);

        double newSpacing[3];
        double newOrigin[3];
        input->GetOrigin(newOrigin);
        input->GetSpacing(newSpacing);

        if (this->zSlice) {
            newdimensions[2] = 1;
            newOrigin[2] += newSpacing[2] * height;
        }
        if (this->xSlice) {
            newdimensions[0] = 1;
            newOrigin[0] += newSpacing[0] * height;
        }
        if (this->ySlice) {
            newdimensions[1] = 1;
            newOrigin[1] += newSpacing[1] * height;
        }

        output->SetDimensions(newdimensions);
        output->SetOrigin(newOrigin);
        output->SetSpacing(newSpacing);
// Iterate all arrays
        for (int arrayIndex = 0; arrayIndex < inPointData->GetNumberOfArrays(); ++arrayIndex) {
            vtkDataArray *inArray = inPointData->GetArray(arrayIndex);

// Allocate an output array
            vtkNew<vtkFloatArray> outArray;
            outArray->SetName(inArray->GetName());

            outArray->SetNumberOfComponents(inArray->GetNumberOfComponents());
            outArray->SetNumberOfTuples(newdimensions[0] * newdimensions[1] * newdimensions[2]);
            outPointData->AddArray(outArray);
// Iterate tuples and scale them
            for (int iz = 0; iz < newdimensions[2]; iz++)
                for (int iy = 0; iy < newdimensions[1]; iy++)
                    for (int ix = 0; ix < newdimensions[0]; ix++) {
                        int tupleIdx = iz * newdimensions[0] * newdimensions[1] + iy * newdimensions[0] + ix;
                        int oldIdx = 0;
                        if (zSlice) {
                            oldIdx = height * dimensions[0] * dimensions[1] + iy * dimensions[0] + ix;
                        }
                        if (xSlice) {
                            oldIdx = iz * dimensions[0] * dimensions[1] + iy * dimensions[0] + height;
                        }
                        if (ySlice) {
                            oldIdx = iz * dimensions[0] * dimensions[1] + height * dimensions[0] + ix;
                        }
                        double *invalue = inArray->GetTuple(oldIdx);
                        outArray->SetTuple(tupleIdx, invalue);
                    }
        }
        if (inPointData->GetScalars())
            outPointData->SetActiveScalars(inPointData->GetScalars()->GetName());
        if (inPointData->GetVectors())
            outPointData->SetActiveVectors(inPointData->GetVectors()->GetName());
        return 1;
    }


private:
    bool zSlice = true;
    bool xSlice = false;
    bool ySlice = false;

    int height = 0;

    BetterVtkSlicer(const BetterVtkSlicer &) = delete;

    void operator=(const BetterVtkSlicer &) = delete;
};

vtkStandardNewMacro(BetterVtkSlicer);

#endif //VIS2021_BETTERVTKSLICER_H
