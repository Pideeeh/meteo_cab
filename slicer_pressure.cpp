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
            height += dataSpacing;

            pressed = true;
        }

        if (key == "Down") {
            height -= dataSpacing;
            pressed = true;
        }
        if (key == "Right") {
            value_lower += 100;
            value_upper += 100;
            pressed = true;
        }

        // Handle a "normal" key
        if (key == "Left") {
            value_lower -= 100;
            value_upper -= 100;
            pressed = true;
        }

        if (pressed) {
            vtkNew<vtkPlane> plane;
            plane->SetOrigin(5, 49, height);
            plane->SetNormal(0, 0, 1);

            sliceMapper->SetSlicePlane(plane);

            origin[2] = height;
            slicer->SetOutputOrigin(origin);
            slicer->Update();
            contourFilter->DebugOn();
            contourFilter->GenerateValues(10, value_lower, value_upper);
            contourFilter->Update();
            contourMapper->Update();
            sliceActor->Update();


            double *pos = contourActor->GetPosition();
            contourActor->SetPosition(pos);

            renderWindow->Render();
        }



        // Forward events
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

    double value_lower = 75000.0;
    double value_upper = 100000.0;
    double height = 0.0;

    vtkSmartPointer<vtkImageReslice> slicer;
    vtkSmartPointer<vtkImageData> data;
    vtkSmartPointer<vtkContourFilter> contourFilter;
    vtkSmartPointer<vtkPolyDataMapper> contourMapper;
    vtkSmartPointer<vtkActor> contourActor;
    vtkSmartPointer<vtkImageActor> sliceActor;
    vtkSmartPointer<vtkImageResliceMapper> sliceMapper;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
};


void CreateColorImage(vtkImageData *input, vtkImageData *image) {
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
    int dimz = image->GetDimensions()[2];

    for (unsigned int x = 0; x < dimx; x++) {
        for (unsigned int y = 0; y < dimy; y++) {
            for (unsigned int z = 0; z < dimz; z++) {
                float *inputColor =
                        static_cast<float *>(input->GetScalarPointer(x, y, z));
                unsigned char *pixel =
                        static_cast<unsigned char *>(image->GetScalarPointer(x, y, z));

                double t = ((inputColor[0] - minVal) / maxVal);
                double color[3];
                ctf->GetColor(t, color);
                for (auto j = 0; j < 3; ++j) {
                    pixel[j] = (unsigned char) (color[j] * 255);
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

    vtkNew<vtkImageData> data_col;
    CreateColorImage(data, data_col);

    vtkNew<vtkImageReslice> reslicer;
    double origin[3];
    data->GetOrigin(origin);

    reslicer->SetOutputOrigin(origin);
    //Todo make this not dependent on a fixed size
    reslicer->SetOutputExtent(0, 199, 0, 199, 0, 0);
    reslicer->SetInputData(data);
    reslicer->Update();


    vtkNew<vtkPlane> plane;
    plane->SetOrigin(5, 49, 19998);
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


    // Create an actor for the contours
    vtkNew<vtkActor> contourActor;
    contourActor->SetMapper(contourMapper);
    contourActor->GetProperty()->SetLineWidth(0.5);

    // Create the outline
    vtkNew<vtkOutlineFilter> outlineFilter;
    outlineFilter->SetInputData(data);

    vtkNew<vtkPolyDataMapper> outlineMapper;
    outlineMapper->SetInputConnection(outlineFilter->GetOutputPort());
    vtkNew<vtkActor> outlineActor;
    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(colors->GetColor3d("Gray").GetData());
    outlineActor->GetProperty()->SetLineWidth(3);

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
    style->slicer = reslicer;
    style->contourFilter = contourFilter;
    style->contourMapper = contourMapper;
    style->contourActor = contourActor;
    style->sliceActor = sliceActor;
    style->sliceMapper = sliceMapper;
    style->height = 19998;
    style->renderWindow = renderWindow;

    style->SetCurrentRenderer(renderer);

    interactor->SetInteractorStyle(style);

    renderWindow->SetSize(500, 500);
    renderWindow->Render();
    interactor->Start();
    return EXIT_SUCCESS;
}