//
// Created by peter on 05.05.2022.
//

#ifndef VIS2021_UTILS_H
#define VIS2021_UTILS_H


#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkNamedColors.h"
#include "vtkColorTransferFunction.h"
#include "vtkXMLImageDataReader.h"

void log(const vtkSmartPointer<vtkObject> &obj){
    obj->PrintSelf(cout,vtkIndent(1));
}



vtkSmartPointer<vtkColorTransferFunction> blueToOrangeTransferFunction(){
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkColorTransferFunction> ctf;
    ctf->SetScaleToLinear();
    ctf->AddRGBPoint(0.0, colors->GetColor3d("MidnightBlue").GetRed(),
                     colors->GetColor3d("MidnightBlue").GetGreen(),
                     colors->GetColor3d("MidnightBlue").GetBlue());
    ctf->AddRGBPoint(0.5, colors->GetColor3d("Gainsboro").GetRed(),
                     colors->GetColor3d("Gainsboro").GetGreen(),
                     colors->GetColor3d("Gainsboro").GetBlue());
    ctf->AddRGBPoint(1.0, colors->GetColor3d("DarkOrange").GetRed(),
                     colors->GetColor3d("DarkOrange").GetGreen(),
                     colors->GetColor3d("DarkOrange").GetBlue());
    return ctf;
}



void
GetSlicewiseColorField(vtkImageData *input, std::vector<float> &localmins, std::vector<float> &localmaxs,
                 vtkImageData *image) {
    image->SetDimensions(input->GetDimensions());
    image->SetSpacing(input->GetSpacing());
    image->SetOrigin(input->GetOrigin());

    auto ctf = blueToOrangeTransferFunction();

    double range[2];
    input->GetScalarRange(range);

    double maxVal = range[1];
    double minVal = range[0];

    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);


    int dimx = image->GetDimensions()[0];
    int dimy = image->GetDimensions()[1];
    const int dimz = image->GetDimensions()[2];

    float mindiff = maxVal;
    float thres = 10000;
    localmins.resize(dimz);
    localmaxs.resize(dimz);


    for (unsigned int z = 0; z < dimz; z++) {
        float max = minVal;
        float min = maxVal;
        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                auto *inputColor =
                        static_cast<float *>(input->GetScalarPointer(x, y, z));
                max = std::max(max, inputColor[0]);
                min = std::min(min, inputColor[0]);
                mindiff = std::min(max - min, mindiff);
            }
        }
        localmins[z] = min;
        localmaxs[z] = max;

        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                float *inputColor =
                        static_cast<float *>(input->GetScalarPointer(x, y, z));
                unsigned char *pixel =
                        static_cast<unsigned char *>(image->GetScalarPointer(x, y, z));

                double t = ((inputColor[0] - min) / std::max((max - min), thres));
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char) (color[j] * 255);
                }
            }
        }
    }
}

void getColorCorrespondingTovalue(double val, double &r, double &g, double &b, double max, double min) {
    double range = max - min;
    static const int numColorNodes = 9;
    double color[numColorNodes][3] =
            {
                    0, 0.7, 0, // green
                    0, 0.6, 0.8,
                    0, 0.4000, 0.6745,
                    0, 0.5000, 0.7745,
                    0, 0.5000, 0.8,
                    0, 0.6000, 0.85,
                    0, 0.6000, 0.9,
                    0, 0.7000, 0.95,
                    0, 0.7000, 1  // Blue
            };

    for (int i = 0; i < (numColorNodes - 1); i++) {
        double currFloor = min + ((double) i / (numColorNodes - 1)) * range;
        double currCeil = min + ((double) (i + 1) / (numColorNodes - 1)) * range;

        if ((val >= currFloor) && (val <= currCeil)) {
            double currFraction = (val - currFloor) / (currCeil - currFloor);
            r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
            g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
            b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
        }
    }
}


void CreateGlobalColorImage(vtkImageData* input, vtkImageData* image) {
    image->SetDimensions(input->GetDimensions());
    image->SetSpacing(input->GetSpacing());
    image->SetOrigin(input->GetOrigin());

    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkColorTransferFunction> ctf;
    ctf->SetScaleToLinear();
    ctf->AddRGBPoint(0.0, colors->GetColor3d("MidnightBlue").GetRed(),
                     colors->GetColor3d("MidnightBlue").GetGreen(),
                     colors->GetColor3d("MidnightBlue").GetBlue());
    ctf->AddRGBPoint(0.5, colors->GetColor3d("Gainsboro").GetRed(),
                     colors->GetColor3d("Gainsboro").GetGreen(),
                     colors->GetColor3d("Gainsboro").GetBlue());
    ctf->AddRGBPoint(1.0, colors->GetColor3d("DarkOrange").GetRed(),
                     colors->GetColor3d("DarkOrange").GetGreen(),
                     colors->GetColor3d("DarkOrange").GetBlue());

    double range[2];
    input->GetScalarRange(range);

    double maxVal = range[1];
    double minVal = range[0];
    double globalmin = minVal;
    double globalmax = maxVal;

    cout << globalmin << endl;
    cout << globalmax << endl;

    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);


    int dimx = image->GetDimensions()[0];
    int dimy = image->GetDimensions()[1];
    const int dimz = image->GetDimensions()[2];

    for (unsigned int z = 0; z < dimz; z++) {
        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                float* inputColor =
                        static_cast<float*>(input->GetScalarPointer(x, y, z));
                unsigned char* pixel =
                        static_cast<unsigned char*>(image->GetScalarPointer(x, y, z));

                double t = ((inputColor[0] - minVal) / (maxVal - minVal));
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char)(color[j] * 255);
                }
            }
        }
    }
}



#endif //VIS2021_UTILS_H
