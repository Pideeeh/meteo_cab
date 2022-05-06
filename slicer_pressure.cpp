#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkColorTransferFunction.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageReslice.h>
#include <vtkPlane.h>
#include <vtkImageActor.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageResliceMapper.h>
#include "dataSetDefinitions.h"
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>

bool horizontal = true;
double globalmin, globalmax;

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera {

public:
    static KeyPressInteractorStyle *New() {
        return new KeyPressInteractorStyle();
    }


vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

public:
    virtual void OnKeyPress() override {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        std::string key = rwi->GetKeySym();
        double dataSpacing = data->GetSpacing()[2];
        double origin[3];
        slicer->GetOutputOrigin(origin);
        bool pressed = false;
        // Handle an arrow key
        if (key == "Up") {
            if (horizontal) {
                height += dataSpacing;
            }
            else {
                vertical_height += 20*dataSpacing;
            }
            pressed = true;
        }

        if (key == "Down") {
            if (horizontal) {
                height -= dataSpacing;
            }
            else {
                vertical_height -= 20*dataSpacing;
            }
            pressed = true;
        }
        if (key == "Right") {
            value_lower += 100;
            value_upper += 100;
            pressed = true;
        }
        if (key == "h") {
            horizontal = !horizontal;
            if (horizontal) {
                cout << "horizontal slice" << endl;
                sliceMapper->SetInputData(data);
            }
            else {
                cout << "vertical slice" << endl;
                sliceMapper->SetInputData(data_global);
            }
        }
        if (key == "Left") {
            value_lower -= 100;
            value_upper -= 100;
            pressed = true;
        }

        if (pressed) {
            vtkNew<vtkPlane> plane;
            plane->SetOrigin(5, vertical_height, height);
            if (horizontal) {
                plane->SetNormal(0, 0, 1);
            }
            else {
                plane->SetNormal(0,1,0);
            }
            

            sliceMapper->SetSlicePlane(plane);

            origin[2] = height;
            slicer->SetOutputOrigin(origin);
            slicer->Update();

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

            // ----------------------------------------------------------------
            // Create a lookup table to share between the mapper and the scalar bar
            // ----------------------------------------------------------------
            float min = localmins[(height-19998)*100];//initialization on values of our height data
            float max = localmaxs[(height - 19998) * 100];
            static const double numColors = 50;
            vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
            lookupTable->SetScaleToLinear();
            lookupTable->SetNumberOfTableValues(numColors);
            double r, g, b;
            for (int i = 0; i < numColors; i++) {
                double val = ((double)i / numColors);
                double color[3];
                ctf->GetColor(val, color);
                lookupTable->SetTableValue(i, color[0], color[1], color[2]);
            }
            lookupTable->Build();
            if (horizontal) {
                lookupTable->SetTableRange(min, max);
            }
            else {
                lookupTable->SetTableRange(globalmin, globalmax);
            }
            // ----------------------------------------------------------------
            // Create a scalar bar actor for the colormap
            // ----------------------------------------------------------------
            legend->SetLookupTable(lookupTable);
            legend->SetNumberOfLabels(3);
            legend->SetTitle("pressure");
            legend->SetVerticalTitleSeparation(6);
            legend->GetPositionCoordinate()->SetValue(0.88, 0.1);
            legend->SetWidth(0.1);

            renderer->AddActor2D(legend);

            contourFilter->DebugOn();
            if (horizontal) {
                contourFilter->GenerateValues(10, min, max);
            }
            else {
                contourFilter->GenerateValues(0, value_lower, value_upper);
            }
            contourFilter->Update();
            contourMapper->Update();
            sliceActor->Update();

            double *pos = contourActor->GetPosition();
            contourActor->SetPosition(pos);
            //cout << contourActor->GetProperty()->GetColor()[0] << contourActor->GetProperty()->GetColor()[1] << contourActor->GetProperty()->GetColor()[1] << endl;

            renderWindow->Render();

        }



        // Forward events
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

    double value_lower = 75000.0;
    double value_upper = 100000.0;
    double height = 0.0;
    double vertical_height = 0.0;

    vtkSmartPointer<vtkImageReslice> slicer;
    vtkSmartPointer<vtkImageData> data, data_global;
    vtkSmartPointer<vtkContourFilter> contourFilter;
    vtkSmartPointer<vtkPolyDataMapper> contourMapper;
    vtkSmartPointer<vtkActor> contourActor;
    vtkSmartPointer<vtkImageActor> sliceActor;
    vtkSmartPointer<vtkImageResliceMapper> sliceMapper;
    vtkSmartPointer<vtkScalarBarActor> legend;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    std::vector<float> localmins;
    std::vector<float> localmaxs;
};


void CreateColorImage(vtkImageData *input, std::vector<float> &localmins, std::vector<float>& localmaxs, vtkImageData *image) {
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

    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);


    int dimx = image->GetDimensions()[0];
    int dimy = image->GetDimensions()[1];
    const int dimz = image->GetDimensions()[2];

    float mindiff = maxVal;
    float thres = 1;
    localmins.resize(dimz);
    localmaxs.resize(dimz);


    for (unsigned int z = 0; z < dimz; z++) {
        float max = minVal;
        float min = maxVal;
        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int x = 0; x < dimx; x++) {
                float* inputColor =
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
                float *inputColor =
                        static_cast<float *>(input->GetScalarPointer(x, y, z));
                unsigned char *pixel =
                        static_cast<unsigned char *>(image->GetScalarPointer(x, y, z));

                double t = ((inputColor[0] - min) / std::max((max-min),thres));
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char) (color[j] * 255);
                }
            }
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
    globalmin = minVal;
    globalmax = maxVal;

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


int main(int, char *[]) {
    vtkNew<vtkNamedColors> colors;

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(getDataPath("/data/press_full.vti").c_str());
    reader->Update();
    vtkSmartPointer<vtkImageData> data = reader->GetOutput();
    data->GetPointData()->SetActiveScalars("pres");

    std::vector<float> localmins;
    std::vector<float> localmaxs;

    vtkNew<vtkImageData> data_col,data_col_global;
    CreateColorImage(data, localmins, localmaxs, data_col);
    CreateGlobalColorImage(data, data_col_global);

    vtkNew<vtkImageReslice> reslicer;
    double origin[3];
    data->GetOrigin(origin);

    reslicer->SetOutputOrigin(origin);
    //Todo make this not dependent on a fixed size
    reslicer->SetOutputExtent(0, 199, 0, 199, 0, 0);
    reslicer->SetInputData(data);
    reslicer->Update();


    vtkNew<vtkPlane> plane;
    plane->SetOrigin(origin);
    plane->SetNormal(0, 0, 1);

    vtkNew<vtkImageResliceMapper> sliceMapper;
    sliceMapper->SetInputData(data_col);
    sliceMapper->SetSlicePlane(plane);


    vtkNew<vtkImageActor> sliceActor;
    sliceActor->SetMapper(sliceMapper);

    // Create an isosurface
    vtkNew<vtkContourFilter> contourFilter;
    contourFilter->SetInputConnection(reslicer->GetOutputPort());
    contourFilter->GenerateValues(10, 75000, 100000); // (numContours, rangeStart, rangeEnd)

// Map the contours to graphical primitives
    vtkNew<vtkPolyDataMapper> contourMapper;
    contourMapper->SetInputConnection(contourFilter->GetOutputPort());
    contourMapper->ScalarVisibilityOn();
    contourMapper->SetScalarRange(100, 1000);

    // Create an actor for the contours
    vtkNew<vtkActor> contourActor;
    contourActor->SetMapper(contourMapper);
    contourActor->GetProperty()->SetLineWidth(1);
    contourActor->GetProperty()->SetColor(colors->GetColor3d("Grey").GetData());

    // Create the outline
    vtkNew<vtkOutlineFilter> outlineFilter;
    outlineFilter->SetInputData(data);

    vtkNew<vtkPolyDataMapper> outlineMapper;
    outlineMapper->SetInputConnection(outlineFilter->GetOutputPort());
    vtkNew<vtkActor> outlineActor;
    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(colors->GetColor3d("Grey").GetData());
    outlineActor->GetProperty()->SetLineWidth(3);


    vtkSmartPointer<vtkScalarBarActor> legend = vtkSmartPointer<vtkScalarBarActor>::New();


    // Visualize
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("IsoContours");

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    renderer->AddActor(contourActor);
    renderer->AddActor(outlineActor);
    renderer->AddActor(sliceActor);
    renderer->SetBackground(0.7, 0.7, 0.7);

    vtkNew<KeyPressInteractorStyle> style;
    style->data = data_col;
    style->data_global = data_col_global;
    style->slicer = reslicer;
    style->contourFilter = contourFilter;
    style->contourMapper = contourMapper;
    style->contourActor = contourActor;
    style->sliceActor = sliceActor;
    style->sliceMapper = sliceMapper;
    style->height = 19998.5;
    style->vertical_height = 48;
    style->renderWindow = renderWindow;
    style->renderer = renderer;
    style->legend = legend;
    style->localmins = localmins;
    style->localmaxs = localmaxs;
    style->SetCurrentRenderer(renderer);

    interactor->SetInteractorStyle(style);

    renderWindow->SetSize(500, 500);
    renderWindow->Render();
    interactor->Start();
    return EXIT_SUCCESS;
}