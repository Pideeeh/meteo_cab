#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCommand.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleUser.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSliderRepresentation.h>
#include <vtkSliderRepresentation3D.h>
#include <vtkSliderWidget.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include "dataSetDefinitions.h"

namespace {
    void CreateData(vtkImageData *data);

    class vtkSliderCallback : public vtkCommand {
    public:
        static vtkSliderCallback *New() {
            return new vtkSliderCallback;
        }

        virtual void Execute(vtkObject *caller, unsigned long, void *) {
            vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget *>(caller);
            double value =
                    static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())
                            ->GetValue();
            this->ContourFilter->GenerateValues(1, value, value);
        }

        vtkSliderCallback() : ContourFilter(NULL) {
        }

        vtkContourFilter *ContourFilter;
    };
} // namespace

int main(int, char *[]) {
    vtkNew<vtkNamedColors> colors;
    //vtkNew<vtkImageData> data;
    //CreateData(data);
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(getDataPath("/data/pressure_slice.vti").c_str());
    reader->Update();
    vtkSmartPointer<vtkImageData> data = reader->GetOutput();

    // Create an isosurface
    vtkNew<vtkContourFilter> contourFilter;
    contourFilter->SetInputData(data);
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
    renderer->SetBackground(1, 1, 1);

    vtkNew<vtkSliderRepresentation3D> sliderRep;
    sliderRep->SetMinimumValue(79000.0);
    sliderRep->SetMaximumValue(100000.0);
    sliderRep->SetValue(80000.0);
    sliderRep->SetTitleText("Contour value");
    sliderRep->SetPoint1InWorldCoordinates(-50, -20, 0);
    sliderRep->SetPoint2InWorldCoordinates(-5, -20, 0);
    sliderRep->SetSliderWidth(0.2);
    sliderRep->SetLabelHeight(0.2);

    vtkNew<vtkSliderWidget> sliderWidget;
    sliderWidget->SetInteractor(interactor);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();

    vtkNew<vtkSliderCallback> callback;
    callback->ContourFilter = contourFilter;

    sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);

    vtkNew<vtkInteractorStyleUser> style;
    interactor->SetInteractorStyle(style);

    renderWindow->SetSize(500, 500);
    renderWindow->Render();
    renderWindow->Render();
    interactor->Start();

    return EXIT_SUCCESS;
}

namespace {
    void CreateData(vtkImageData *data) {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(
                "C:/Users/annik/Documents/fs22/scientific visualization/vis22_prog/meteo_cab/data/pres_final.vti");
        reader->Update();
        vtkSmartPointer<vtkImageData> dat = reader->GetOutput();
        int *dims = dat->GetDimensions();
        vtkSmartPointer<vtkPointData> pts = dat->GetPointData();
        vtkSmartPointer<vtkAbstractArray> arr = pts->GetAbstractArray("pres");
        vtkSmartPointer<vtkFloatArray> float_array = vtkFloatArray::SafeDownCast(arr);
        data->SetExtent(0, 255, 0, 255, 0, 0);
        int *extent = data->GetExtent();
        data->AllocateScalars(VTK_DOUBLE, 1);
        for (int y = extent[2]; y <= extent[3]; y++) {
            for (int x = extent[0]; x <= extent[1]; x++) {
                double *pixel = static_cast<double *>(data->GetScalarPointer(x, y, 0));
                int idx = y * dims[0] + x;
                auto v = static_cast<double *>(float_array->GetTuple(idx));
                pixel[0] = v[0];
            }
        }

        vtkNew<vtkXMLImageDataWriter> writer;
        writer->SetFileName("data.vti");
        writer->SetInputData(data);
        writer->Write();
    }
} // namespace
