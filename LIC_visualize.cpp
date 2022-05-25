#include <random>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include "dataSetDefinitions.h"
#include "vtkImageCast.h"
#include "vtkImageMapper.h"
#include "vtkActor2D.h"
#include <vtkCharArray.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPiecewiseFunction.h>
#include "BetterVtkSlicer.h"

#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkContourFilter.h>
#include <vtkPNGReader.h>
#include <vtkInteractorStyleUser.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>

namespace {
    void CreateData(vtkImageData* data);

    class vtkSliderCallback : public vtkCommand {
    public:
        static vtkSliderCallback* New() {
            return new vtkSliderCallback;
        }

        virtual void Execute(vtkObject* caller, unsigned long, void*) {
            vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
            double value =
                static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())
                ->GetValue();

            // Load the image
            vtkNew<vtkPNGReader> reader;
            int idx = int(((value - 183) / (20000 - 183)) * 199);
            reader->SetFileName(getDataPath("/data/lic" + std::to_string(idx) + ".png").c_str());
            reader->Update();

            // Give image to the mapper
            img = reader->GetOutput();
            img_map->SetInputData(img);
        }

        vtkSliderCallback() : img_map(NULL), img(NULL) {
        }

        vtkImageMapper *img_map;
        vtkImageData* img;
    };
} // namespace

int main(int, char* []) {
    vtkNew<vtkNamedColors> colors;
    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
    vtkNew<vtkImageMapper> img_map;
    vtkNew<vtkActor2D> image_actor;
    image_actor->SetMapper(img_map);
    img_map->SetColorWindow(255);
    img_map->SetColorLevel(127.5);

    // Initialization of img and img_map for zero level values
    vtkNew<vtkPNGReader> reader;;
    reader->SetFileName(getDataPath("/data/lic" + std::to_string(0) + ".png").c_str());
    reader->Update();
    img = reader->GetOutput();
    img_map->SetInputData(img);

    // Visualize
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("LIC texture");

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(renderWindow);

    renderer->AddActor(image_actor);
    renderer->SetBackground(1, 1, 1);

    vtkNew<vtkSliderRepresentation2D> sliderRep;
    sliderRep->SetMinimumValue(183);
    sliderRep->SetMaximumValue(20000);
    sliderRep->SetValue(183);
    sliderRep->SetValue(0);
    sliderRep->SetTitleText("Altitude (meters)");

    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(160, 550);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(560, 550);
    sliderRep->SetSliderWidth(0.02);
    sliderRep->SetLabelHeight(0.02);
    sliderRep->SetSliderLength(0.03);
    sliderRep->GetTubeProperty()->SetColor(
        colors->GetColor3d("Red").GetData());
    sliderRep->GetSliderProperty()->SetColor(
        colors->GetColor3d("Blue").GetData());
    sliderRep->GetLabelProperty()->SetColor(
        colors->GetColor3d("Black").GetData());
    sliderRep->GetTitleProperty()->SetColor(
        colors->GetColor3d("Black").GetData());

    vtkNew<vtkSliderWidget> sliderWidget;
    sliderWidget->SetInteractor(interactor);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();

    vtkNew<vtkSliderCallback> callback;
    callback->img_map = img_map;
    callback->img = img;

    sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);

    vtkNew<vtkInteractorStyleUser> style;
    interactor->SetInteractorStyle(style);

    renderWindow->SetSize(720, 600);
    renderWindow->Render();
    renderWindow->Render();
    interactor->Start();

    return EXIT_SUCCESS;
}