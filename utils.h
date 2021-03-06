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

void log(const vtkSmartPointer<vtkObject>& obj) {
    obj->PrintSelf(cout, vtkIndent(1));
}


vtkSmartPointer<vtkColorTransferFunction> blueToOrangeTransferFunction() {
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
GetSlicewiseColorField(vtkImageData* input, std::vector<float>& localmins, std::vector<float>& localmaxs,
    vtkImageData* image) {
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
                auto* inputColor =
                    static_cast<float*>(input->GetScalarPointer(x, y, z));
                max = std::max(max, inputColor[0]);
                min = std::min(min, inputColor[0]);
                mindiff = std::min(max - min, mindiff);
            }
        }
        localmins[z] = min;
        localmaxs[z] = max;

        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                float* inputColor =
                    static_cast<float*>(input->GetScalarPointer(x, y, z));
                unsigned char* pixel =
                    static_cast<unsigned char*>(image->GetScalarPointer(x, y, z));

                double t = ((inputColor[0] - min) / std::max((max - min), thres));
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char)(color[j] * 255);
                }
            }
        }
    }
}

void getColorCorrespondingTovalue(double val, double& r, double& g, double& b, double max, double min) {
    double range = max - min;
    static const int numColorNodes = 7;
    double color[numColorNodes][3] =
    {
            103. / 255, 130. / 255, 113. / 255,//green
            137. / 255, 166. / 255, 148. / 255,//a bit lighter green
            185. / 255, 180. / 255, 150. / 255,//brown
            185. / 255 + 10. / 255, 180. / 255 + 10. / 255, 150. / 255 + 10. / 255,//lighter brown
            185. / 255 + 30. / 255, 180. / 255 + 30. / 255, 150. / 255 + 30. / 255,//lighter brown
            1, 1, 1,//white
            1, 1, 1,//white
    };

    for (int i = 0; i < (numColorNodes - 1); i++) {
        double currFloor = min + ((double)i / (numColorNodes - 1)) * range;
        double currCeil = min + ((double)(i + 1) / (numColorNodes - 1)) * range;

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

    double range[2];
    input->GetScalarRange(range);

    ctf->SetScaleToLinear();
    /*ctf->AddRGBPoint(, colors->GetColor3d("MidnightBlue").GetRed(),
                     colors->GetColor3d("MidnightBlue").GetGreen(),
                     colors->GetColor3d("MidnightBlue").GetBlue());*/
    ctf->AddRGBPoint(0.0, colors->GetColor3d("Gainsboro").GetRed(),
        colors->GetColor3d("Gainsboro").GetGreen(),
        colors->GetColor3d("Gainsboro").GetBlue());
    ctf->AddRGBPoint(100, colors->GetColor3d("DarkOrange").GetRed(),
        colors->GetColor3d("DarkOrange").GetGreen(),
        colors->GetColor3d("DarkOrange").GetBlue());

    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    int dimx = image->GetDimensions()[0];
    int dimy = image->GetDimensions()[1];
    const int dimz = image->GetDimensions()[2];

    for (unsigned int z = 0; z < dimz; z++) {
        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                double* inputColor =
                    static_cast<double*>(input->GetScalarPointer(x, y, z));
                unsigned char* pixel =
                    static_cast<unsigned char*>(image->GetScalarPointer(x, y, z));

                double t = abs(inputColor[0]);
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char)(color[j] * 255);
                }
            }
        }
    }
}

vtkSmartPointer<vtkPiecewiseFunction>
GetDivergenceOpacityFunction(double minValue, double maxValue, double alpha, double transparentFrom,
    double transparentTo) {
    vtkSmartPointer<vtkPiecewiseFunction> opacity = vtkNew<vtkPiecewiseFunction>();
    opacity->AddPoint(-100, alpha);
    opacity->AddPoint(100, alpha);
    opacity->AddPoint(transparentFrom, alpha);
    opacity->AddPoint(transparentFrom + 0.01, 0);
    opacity->AddPoint(transparentTo, alpha);
    opacity->AddPoint(transparentTo - 0.01, 0);
    return opacity;
}


vtkSmartPointer<vtkColorTransferFunction>
GetDivergenceColorFunction(double minValue, double maxValue) {
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkNew<vtkColorTransferFunction>();

    ctf->AddRGBPoint(minValue, 1, 0, 0);
    ctf->AddRGBPoint(-1, 1, 0, 0);

    ctf->AddRGBPoint(1, 0, 0, 1);
    ctf->AddRGBPoint(maxValue, 0, 0, 1);

    return ctf;
}

vtkSmartPointer<vtkPiecewiseFunction>
GetSimpleOpacityFunction(double maxValue, double alpha, double transparentTo) {
    vtkSmartPointer<vtkPiecewiseFunction> opacity = vtkNew<vtkPiecewiseFunction>();
    //opacity->AddPoint(minValue, alpha);
    opacity->AddPoint(maxValue, alpha);
    //opacity->AddPoint(transparentFrom, alpha);
    //opacity->AddPoint(transparentFrom + 0.01, 0);
    opacity->AddPoint(transparentTo, alpha);
    opacity->AddPoint(transparentTo - 0.00001, 0);
    opacity->AddPoint(transparentTo - 1, 0);

    return opacity;
}

vtkSmartPointer<vtkColorTransferFunction>
GetSimpleColorFunction(double maxValue, double r, double g, double b) {
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkNew<vtkColorTransferFunction>();

    ctf->AddRGBPoint(0, r, g, b);
    //ctf->AddRGBPoint(-1, 1, 0, 0);
    //ctf->AddRGBPoint(1, 0, 0, 1);
    ctf->AddRGBPoint(maxValue, r, g, b);

    return ctf;
}
#endif //VIS2021_UTILS_H
