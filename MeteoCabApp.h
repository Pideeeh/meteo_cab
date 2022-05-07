//
// Created by peter on 05.05.2022.
//

#ifndef VIS2021_METEOCABAPP_H
#define VIS2021_METEOCABAPP_H

#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkNamedColors.h"
#include "vtkColorTransferFunction.h"
#include "vtkXMLImageDataReader.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkScalarsToColors.h"
#include "vtkPlane.h"
#include "vtkImageActor.h"
#include "vtkImageResliceMapper.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkScalarBarActor.h"
#include "vtkOutlineFilter.h"
#include "vtkProperty.h"
#include "vtkOBJReader.h"
#include "vtkPlaneWidget.h"
#include "vtkFloatArray.h"

#include "utils.h"
#include "dataSetDefinitions.h"
#include "PlaneWidgetInteraction.h"

using namespace std;



class MeteoCabApp : public vtkInteractorStyleTrackballCamera {
    static MeteoCabApp *New() {
        return new MeteoCabApp();
    }

vtkTypeMacro(MeteoCabApp, vtkInteractorStyleTrackballCamera);
    vtkNew<vtkNamedColors> colors;

/*
 * Pressure Variables
 * */
    bool pressureVisible = true;
    vtkSmartPointer<vtkImageData> pressureData;
    vtkSmartPointer<vtkImageData> pressureDataColors;
    vtkSmartPointer<vtkImageActor> pressureSliceActor;
    vtkSmartPointer<vtkActor> contourActor;
    vtkSmartPointer<vtkActor> outlineActor;
    vtkSmartPointer<vtkContourFilter> pressureContourFilter;
    vtkSmartPointer<vtkPolyDataMapper> pressureContourMapper;
    vtkSmartPointer<vtkImageResliceMapper> pressureSliceMapper;

    vtkSmartPointer<vtkScalarBarActor> pressureLegend;

    vtkSmartPointer<vtkImageReslice> pressureReslicer;
    vtkSmartPointer<PlaneWidgetInteraction> planeWidgetInteraction;

    std::vector<float> pressure_slice_minima;
    std::vector<float> pressure_slice_maxima;

/*
 * Heightmap Variables
 * */
    bool heightmapVisible = false;
    vtkSmartPointer<vtkActor> heightmapActor;
    vtkSmartPointer<vtkScalarBarActor> heightmapLegend;

/*
 * Rendering Variables
 * */
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

/*
 * UI Variables
 * */
    vtkSmartPointer<vtkPlaneWidget> planeWidget;
public:
    void StartVisualization() {
        renderer = vtkNew<vtkRenderer>();
        renderWindow = vtkNew<vtkRenderWindow>();
        renderWindow->AddRenderer(renderer);
        renderWindow->SetWindowName("MeteoCab");

        interactor = vtkNew<vtkRenderWindowInteractor>();
        interactor->SetRenderWindow(renderWindow);

        InitializePlaneWidget();

        renderer->AddActor(contourActor);
        renderer->AddActor(outlineActor);
        renderer->AddActor(pressureSliceActor);
        renderer->AddActor(heightmapActor);
        renderer->SetBackground(0.7, 0.7, 0.7);
        interactor->SetInteractorStyle(this);
        renderWindow->SetSize(500, 500);
        renderWindow->Render();
        interactor->Start();
    }

    void InitializePressureSlicer() {
        vtkNew<vtkXMLImageDataReader> reader;
        reader->SetFileName(getDataPath("/data/press_full.vti").c_str());
        reader->Update();
        this->pressureData = reader->GetOutput();
        this->pressureData->GetPointData()->SetActiveScalars("pres");
        pressureDataColors = vtkNew<vtkImageData>();

        GetSlicewiseColorField(pressureData, pressure_slice_minima, pressure_slice_maxima, pressureDataColors);

        pressureReslicer = vtkNew<vtkImageReslice>();
        double origin[3];
        pressureData->GetOrigin(origin);

        pressureReslicer->SetOutputOrigin(origin);
        //Todo make this not dependent on a fixed size
        int extent[6];
        pressureData->GetExtent(extent);
        extent[5] = 0;
        pressureReslicer->SetOutputExtent(extent);
        pressureReslicer->SetInputData(pressureData);
        pressureReslicer->Update();

        vtkNew<vtkPlane> plane;
        plane->SetOrigin(5, 49, 19998);
        plane->SetNormal(0, 0, 1);

        pressureSliceMapper = vtkNew<vtkImageResliceMapper>();
        pressureSliceMapper->SetInputData(pressureDataColors);
        pressureSliceMapper->SetSlicePlane(plane);

        pressureSliceActor = vtkNew<vtkImageActor>();
        pressureSliceActor->SetMapper(pressureSliceMapper);

        // Create an isosurface
        pressureContourFilter = vtkNew<vtkContourFilter>();
        pressureContourFilter->SetInputConnection(pressureReslicer->GetOutputPort());
        pressureContourFilter->GenerateValues(10, 75000, 100000); // (numContours, rangeStart, rangeEnd)

        // Map the contours to graphical primitives
        pressureContourMapper = vtkNew<vtkPolyDataMapper>();
        pressureContourMapper->SetInputConnection(pressureContourFilter->GetOutputPort());


        // Create an actor for the contours
        contourActor = vtkNew<vtkActor>();
        contourActor->SetMapper(pressureContourMapper);
        contourActor->GetProperty()->SetLineWidth(0.5);

        vtkNew<vtkOutlineFilter> outlineFilter;
        outlineFilter->SetInputData(pressureData);

        vtkNew<vtkPolyDataMapper> outlineMapper;
        outlineMapper->SetInputConnection(outlineFilter->GetOutputPort());
        outlineActor = vtkNew<vtkActor>();
        outlineActor->SetMapper(outlineMapper);
        outlineActor->GetProperty()->SetColor(colors->GetColor3d("Grey").GetData());
        outlineActor->GetProperty()->SetLineWidth(3);

        pressureLegend = vtkNew<vtkScalarBarActor>();

        TogglePressure();
        TogglePressure();
    }

    void InitializeHeightmap() {
        // ----------------------------------------------------------------
        // color initialization based on height
        // ----------------------------------------------------------------

        double height_min = 0;//initialization on values of our height data
        double height_max = 0.244379;
        double numColors = 20;
        int elementsx = 1429;//how many x values on grid
        int elementsy = 1556;//how many y
        vtkNew<vtkFloatArray> height_colors;
        height_colors->SetNumberOfValues(elementsx * elementsy);
        // ----------------------------------------------------------------
        // read in text file with height data to decide on color
        // ----------------------------------------------------------------
        std::ifstream myfile(getDataPath("/data/height_values.txt"));
        if (!myfile.is_open()) {
            std::cout << "error opening file" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (unsigned int ix = 0; ix < elementsx; ix++) {
            for (unsigned int iy = 0; iy < elementsy; iy++) {
                std::string str;
                std::getline(myfile, str);
                double h = std::stod(str) * 0.0001;
                // Height value color-coded
                height_colors->SetValue(ix * elementsy + iy, h);
            }
        }
        // ----------------------------------------------------------------
        // Create a lookup table to share between the mapper and the scalar bar
        // ----------------------------------------------------------------
        vtkNew<vtkLookupTable> lookupTable;
        lookupTable->SetScaleToLinear();
        lookupTable->SetNumberOfTableValues(numColors);
        double r, g, b;
        for (int i = 0; i < numColors; i++) {
            double val = height_min + ((double) i / numColors) * (height_max - height_min);
            getColorCorrespondingTovalue(val, r, g, b, height_max, height_min);
            lookupTable->SetTableValue(i, r, g, b);
        }
        lookupTable->Build();
        // ----------------------------------------------------------------
        // Create a scalar bar actor for the colormap
        // ----------------------------------------------------------------
        heightmapLegend = vtkNew<vtkScalarBarActor>();
        heightmapLegend->SetLookupTable(lookupTable);
        heightmapLegend->SetNumberOfLabels(3);
        heightmapLegend->SetTitle("height");
        heightmapLegend->SetVerticalTitleSeparation(6);
        heightmapLegend->GetPositionCoordinate()->SetValue(0.88, 0.1);
        heightmapLegend->SetWidth(0.1);
        // ----------------------------------------------------------------
        // Create a triangulated mapper for mesh read from file, color it and link it to the lookup table
        // ----------------------------------------------------------------
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName(getDataPath("/data/heightfield.obj").c_str());
        reader->Update();//read .obj file
        vtkNew<vtkPolyDataMapper> triangulatedMapper;
        triangulatedMapper->SetInputConnection(reader->GetOutputPort());//take mesh as input
        triangulatedMapper->GetInput()->GetPointData()->SetScalars(height_colors);//color it based on our computations
        triangulatedMapper->SetLookupTable(lookupTable);//link lookup table
        triangulatedMapper->SetScalarRange(height_min, height_max);
        // ----------------------------------------------------------------
        // Create an actor
        // ----------------------------------------------------------------
        heightmapActor = vtkNew<vtkActor>();
        heightmapActor->SetMapper(triangulatedMapper);
        heightmapActor->GetProperty()->SetEdgeVisibility(false);
    }

    void InitializePlaneWidget() {
        planeWidget = vtkNew<vtkPlaneWidget>();
        planeWidget->SetInteractor(interactor);
        planeWidget->SetRepresentationToOutline();

        planeWidgetInteraction = vtkNew<PlaneWidgetInteraction>();
        planeWidgetInteraction->dataSpace = pressureData;
        planeWidgetInteraction->app = this;
        planeWidget->AddObserver(vtkCommand::InteractionEvent, planeWidgetInteraction);
        planeWidget->AddObserver(vtkCommand::EndInteractionEvent, planeWidgetInteraction);

        double center[3];
        pressureData->GetCenter(center);
        double origin[3];
        pressureData->GetOrigin(origin);
        center[2] = origin[2];
        double bounds[6];
        pressureData->GetBounds(bounds);
        planeWidget->SetCenter(center);
        planeWidget->SetOrigin(bounds[0], bounds[2], bounds[4]);
        planeWidget->SetPoint1(bounds[0], bounds[3], bounds[4]);
        planeWidget->SetPoint2(bounds[1], bounds[2], bounds[4]);


        //planeWidget->PlaceWidget(bounds);
        planeWidget->On();
    }

    void Launch() {
        InitializePressureSlicer();
        InitializeHeightmap();

        StartVisualization();
    }

    void UpdatePressureSlicer() {

        int dimensions[3];
        pressureData->GetDimensions(dimensions);

        double dataOrigin[3];
        pressureData->GetOrigin(dataOrigin);

        double dataSpacing = pressureData->GetSpacing()[2];
        double origin[3];
        pressureReslicer->GetOutputOrigin(origin);

        vtkNew<vtkPlane> plane;
        planeWidget->GetPlane(plane);

        double planeOrigin[3];
        plane->GetOrigin(planeOrigin);

        double sliceHeight = planeOrigin[2];
        double current_slice = (sliceHeight - dataOrigin[2]) / dataSpacing;

        pressureSliceMapper->SetSlicePlane(plane);

        origin[2] = sliceHeight;
        pressureReslicer->SetOutputOrigin(origin);
        pressureReslicer->Update();

        vtkSmartPointer<vtkColorTransferFunction> ctf = blueToOrangeTransferFunction();

        // ----------------------------------------------------------------
        // Create a lookup table to share between the mapper and the scalar bar
        // ----------------------------------------------------------------
        float min = pressure_slice_minima[current_slice];
        float max = pressure_slice_maxima[current_slice];
        static const double numColors = 10;
        vtkSmartPointer<vtkLookupTable> lookupTable = vtkNew<vtkLookupTable>();
        lookupTable->SetScaleToLinear();
        lookupTable->SetNumberOfTableValues(numColors);
        for (int i = 0; i < numColors; i++) {
            double val = ((double) i / numColors);
            double color[3];
            ctf->GetColor(val, color);
            lookupTable->SetTableValue(i, color[0], color[1], color[2]);
        }
        lookupTable->Build();
        lookupTable->SetTableRange(min, max);
        // ----------------------------------------------------------------
        // Create a scalar bar actor for the colormap
        // ----------------------------------------------------------------
        pressureLegend->SetLookupTable(lookupTable);
        pressureLegend->SetNumberOfLabels(3);
        pressureLegend->SetTitle("pressure");
        pressureLegend->SetVerticalTitleSeparation(6);
        pressureLegend->GetPositionCoordinate()->SetValue(0.88, 0.1);
        pressureLegend->SetWidth(0.1);

        renderer->AddActor2D(pressureLegend);

        pressureContourFilter->GenerateValues(10, pressure_value_lower, pressure_value_upper);
        pressureContourFilter->Update();
        pressureContourMapper->Update();
        pressureSliceActor->Update();

        double *pos = contourActor->GetPosition();
        contourActor->SetPosition(pos);

    }

    void PlaneWidgetUpdated() {
        UpdatePressureSlicer();

        renderWindow->Render();
    }


    void TogglePressure() {
        pressureVisible = !pressureVisible;
        pressureSliceActor->SetVisibility(pressureVisible);
        contourActor->SetVisibility(pressureVisible);
        pressureLegend->SetVisibility(pressureVisible);
    }

    void ToggleHeightmap() {
        heightmapVisible = !heightmapVisible;
        heightmapActor->SetVisibility(heightmapVisible);
        heightmapLegend->SetVisibility(heightmapVisible);

    }

    void HandlePressureInputs(const string &key) {
        // Handle an arrow key
        bool pressed = false;
        if (key == "Right") {
            pressure_value_lower += 100;
            pressure_value_upper += 100;
            pressed = true;
        }

        // Handle a "normal" key
        if (key == "Left") {
            pressure_value_lower -= 100;
            pressure_value_upper -= 100;
            pressed = true;
        }

        if (pressed) {
            UpdatePressureSlicer();
        }
    }

    void OnLeftButtonUp() override {
        if (planeWidgetInteraction->updated) {
            PlaneWidgetUpdated();
            planeWidgetInteraction->updated = false;
        }
        vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
    }

    void OnKeyPress() override {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        string key = rwi->GetKeySym();

        cout << "Key " << key << endl;

        if (key == "1") {
            cout << "Toggle Pressure" << endl;
            TogglePressure();
        }

        if (key == "2") {
            ToggleHeightmap();
        }

        if (pressureVisible) {
            HandlePressureInputs(key);
        }

        renderWindow->Render();


        // Forward events
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

    double pressure_value_lower = 75000.0;
    double pressure_value_upper = 100000.0;


};


#endif //VIS2021_METEOCABAPP_H
