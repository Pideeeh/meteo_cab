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
#include "vtkTransform.h"
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
#include "vtkMetaImageReader.h"
#include "vtkVolumeProperty.h"
#include "vtkEasyTransfer.hpp"
#include "vtkTransformPolyDataFilter.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkMaskPoints.h"
#include "vtkGlyphSource2D.h"
#include "vtkGradientFilter.h"
#include "vtkGlyph2D.h"

#include "utils.h"
#include "dataSetDefinitions.h"
#include "PlaneWidgetInteraction.h"
#include "BetterVtkSlicer.h"

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
    vtkSmartPointer<vtkImageData> slice;
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
    vtkSmartPointer<vtkPolyDataMapper> heightmapDataMapper;

    vtkSmartPointer<vtkLookupTable> heightmapLookupTable;
    vtkSmartPointer<vtkFloatArray> heightmapColors;


/*
 * 2D Wind Variables
 * */
    int windMode = 1;
    bool divergenceVisible = false;
    vtkSmartPointer<vtkImageData> horizontalWindDivergence;
    vtkSmartPointer<vtkImageData> horizontalWindDivergenceColors;
    vtkSmartPointer<vtkImageActor> horizontalWindDivergenceActor;
    vtkSmartPointer<vtkImageResliceMapper> horizontalWindDivergenceSliceMapper;

    vtkSmartPointer<vtkActor> horizontalWindVectorActor;
    vtkSmartPointer<vtkImageData> horizontalWindData;
    vtkSmartPointer<BetterVtkSlicer> horizontalWindSlicer;
    vtkSmartPointer<vtkPolyDataMapper> windMapper;
    vtkSmartPointer<vtkGlyph2D> windGlyphs;
    vtkSmartPointer<vtkGlyphSource2D> windGlyphSource;
    vtkSmartPointer<vtkMaskPoints> windMask;
/*
 * Rendering Variables
 * */
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

    /*
     * Cloud variables
     * */
    bool cloudVisible = false;
    vtkSmartPointer<vtkImageData> cloudData;
    vtkSmartPointer<vtkVolume> cloudActor;
    vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> cloudRaycastMapper;
    vtkSmartPointer<vtkEasyTransfer> cloudTransferFunction;
    vtkSmartPointer<vtkRenderer> cloudTransferRenderer;
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
        renderer->AddActor(cloudActor);
        renderer->AddActor2D(heightmapLegend);
        renderer->AddActor(horizontalWindVectorActor);
        renderer->AddActor(horizontalWindDivergenceActor);
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
        double baseHeight = pressureData->GetOrigin()[2];
        double height_min = baseHeight + 0;//initialization on values of our height data
        double height_max = baseHeight + 0.244379;
        double numColors = 20;
        
        // ----------------------------------------------------------------
        // Create a triangulated mapper for mesh read from file, color it and link it to the lookup table
        // ----------------------------------------------------------------
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName(getDataPath("/data/downsampled_terrain.obj").c_str());
        reader->Update();//read .obj file
        heightmapDataMapper = vtkNew<vtkPolyDataMapper>();
        heightmapDataMapper->SetInputConnection(reader->GetOutputPort());

        heightmapColors = vtkNew<vtkFloatArray>();
        heightmapColors->SetNumberOfValues(229972);
        for (int i = 0; i < 229972; i++) {//go over vertices
            heightmapColors->SetValue(i, heightmapDataMapper->GetInput()->GetPoint(i)[2]+baseHeight);
        }
        // ----------------------------------------------------------------
        // Create a lookup table to share between the mapper and the scalar bar
        // ----------------------------------------------------------------
        heightmapLookupTable = vtkNew<vtkLookupTable>();
        heightmapLookupTable->SetScaleToLinear();
        heightmapLookupTable->SetNumberOfTableValues(numColors);
        double r, g, b;
        for (int i = 0; i < numColors; i++) {
            double val = height_min + ((double)i / numColors) * (height_max - height_min);
            getColorCorrespondingTovalue(val, r, g, b, height_max, height_min);
            heightmapLookupTable->SetTableValue(i, r, g, b);
        }
        heightmapLookupTable->Build();
        // ----------------------------------------------------------------
        // Create a scalar bar actor for the colormap
        // ----------------------------------------------------------------
        heightmapLegend = vtkNew<vtkScalarBarActor>();
        heightmapLegend->SetLookupTable(heightmapLookupTable);
        heightmapLegend->SetNumberOfLabels(3);
        heightmapLegend->SetTitle("height");
        heightmapLegend->SetVerticalTitleSeparation(6);
        heightmapLegend->GetPositionCoordinate()->SetValue(0.85, 0.1);
        heightmapLegend->SetWidth(0.05);

        heightmapDataMapper->GetInput()->GetPointData()->SetScalars(
                heightmapColors);//color it based on our computations
        heightmapDataMapper->SetLookupTable(heightmapLookupTable);//link lookup table
        heightmapDataMapper->SetScalarRange(height_min, height_max);
        // ----------------------------------------------------------------
        // Create an actor
        // ----------------------------------------------------------------
        heightmapActor = vtkNew<vtkActor>();
        heightmapActor->SetMapper(heightmapDataMapper);
        heightmapActor->GetProperty()->SetEdgeVisibility(false);
        double origin[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = pressureData->GetOrigin()[2];
        heightmapActor->SetPosition(origin);
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

        planeWidget->On();
    }

    void InitializeClouds() {

        vtkNew<vtkMetaImageReader> reader;
        reader->SetFileName(getDataPath("/data/cloudyness.mhd").c_str());
        reader->Update();
        cloudData = reader->GetOutput();

        cloudRaycastMapper = vtkNew<vtkOpenGLGPUVolumeRayCastMapper>();
        cloudRaycastMapper->SetInputData(cloudData);

        // create transfer function
        cloudTransferFunction = vtkNew<vtkEasyTransfer>();
        cloudTransferFunction->SetColorUniform(1.0, 1.0, 1.0);        // set initial color map
        cloudTransferFunction->SetColorRange(0, 0.00507);    // set the value range that is mapped to color
        cloudTransferFunction->SetCloudOpacityRange(0, 0.00507, 3, 1000,
                                                    0.9);    // set the value range that is mapped to opacity
        cloudTransferFunction->RefreshImage();

        // assign transfer function to volume properties
        vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        volumeProperty->SetColor(cloudTransferFunction->GetColorTransferFunction());
        volumeProperty->SetScalarOpacity(cloudTransferFunction->GetOpacityTransferFunction());

        // create volume actor and assign mapper and properties
        cloudActor = vtkNew<vtkVolume>();
        cloudActor->SetMapper(cloudRaycastMapper);
        cloudActor->SetProperty(volumeProperty);
        cloudActor->SetVisibility(cloudVisible);

        // get renderer for the transfer function editor or a white background
        cloudTransferRenderer = cloudTransferFunction->GetRenderer();
        cloudTransferRenderer->SetViewport(0.66666, 0, 1, 1);

    }

    void InitializeHorizontalWind() {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/2d_wind_full.vti").c_str());
        reader->Update();
        horizontalWindData = reader->GetOutput();
        horizontalWindData->GetPointData()->SetActiveVectors("2d_velocity");

        double origin[3];
        horizontalWindData->GetOrigin(origin);

        horizontalWindSlicer = vtkNew<BetterVtkSlicer>();
        horizontalWindSlicer->SetZSlice();
        horizontalWindSlicer->SetHeight(0);
        horizontalWindSlicer->SetInputData(horizontalWindData);
        horizontalWindSlicer->Update();

        windMask = vtkNew<vtkMaskPoints>();

        windMask->SetInputData(horizontalWindSlicer->GetOutput());
        windMask->RandomModeOn();
        windMask->SetRandomModeType(2);
        windMask->SetMaximumNumberOfPoints(1000);
        windMask->SetOnRatio(50);
        windMask->Update();

        windGlyphSource = vtkNew<vtkGlyphSource2D>();
        windGlyphSource->SetGlyphTypeToArrow();
        windGlyphSource->FilledOff();

        windGlyphs = vtkNew<vtkGlyph2D>();
        windGlyphs->SetInputConnection(windMask->GetOutputPort());
        windGlyphs->SetSourceConnection(windGlyphSource->GetOutputPort());
        windGlyphs->OrientOn();
        windGlyphs->SetScaleModeToScaleByVector();
        windGlyphs->SetScaleFactor(0.04);
        windGlyphs->Update();

        windMapper = vtkNew<vtkPolyDataMapper>();
        windMapper->SetInputConnection(windGlyphs->GetOutputPort());
        windMapper->Update();

        horizontalWindVectorActor = vtkNew<vtkActor>();
        horizontalWindVectorActor->SetMapper(windMapper);
        horizontalWindVectorActor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
        double position[3];
        horizontalWindVectorActor->GetPosition(position);
        position[2] = horizontalWindData->GetBounds()[4];
        horizontalWindVectorActor->SetPosition(position);
    }

    void ComputeWindDivergence() {

        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/2d_wind_full.vti").c_str());
        reader->Update();

        horizontalWindDivergence = reader->GetOutput();

        vtkNew<vtkGradientFilter> gradientFilter;
        gradientFilter->SetInputConnection(0, reader->GetOutputPort());
        gradientFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "2d_velocity");
        gradientFilter->ComputeDivergenceOn();
        gradientFilter->Update();

        auto div = gradientFilter->GetOutput();
        vtkImageData *input = vtkImageData::SafeDownCast(div);
        input->GetPointData()->SetActiveScalars("Divergence");
        horizontalWindDivergence = input;

        horizontalWindDivergenceColors = vtkNew<vtkImageData>();
        CreateGlobalColorImage(horizontalWindDivergence, horizontalWindDivergenceColors);

        vtkNew<vtkPlane> plane;
        plane->SetOrigin(5, 49, 19998);
        plane->SetNormal(0, 0, 1);

        horizontalWindDivergenceSliceMapper = vtkNew<vtkImageResliceMapper>();
        horizontalWindDivergenceSliceMapper->SetInputData(horizontalWindDivergenceColors);
        horizontalWindDivergenceSliceMapper->SetSlicePlane(plane);

        horizontalWindDivergenceActor = vtkNew<vtkImageActor>();
        horizontalWindDivergenceActor->SetMapper(horizontalWindDivergenceSliceMapper);
    }

    void Launch() {
        InitializePressureSlicer();
        InitializeHeightmap();

        InitializeHorizontalWind();
        ComputeWindDivergence();
        InitializeClouds();
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
        pressureLegend->GetPositionCoordinate()->SetValue(0.93, 0.1);
        pressureLegend->SetWidth(0.05);

        renderer->AddActor2D(pressureLegend);

        pressureContourFilter->GenerateValues(10, pressure_value_lower, pressure_value_upper);
        pressureContourFilter->Update();
        pressureContourMapper->Update();
        pressureSliceActor->Update();

        double *pos = contourActor->GetPosition();
        contourActor->SetPosition(pos);

    }

    void UpdateWindSlice() {
        int dimensions[3];
        horizontalWindData->GetDimensions(dimensions);

        double dataOrigin[3];
        horizontalWindData->GetOrigin(dataOrigin);

        double dataSpacing = horizontalWindData->GetSpacing()[2];
        double origin[3];
        horizontalWindSlicer->GetOutput()->GetOrigin(origin);

        vtkNew<vtkPlane> plane;
        planeWidget->GetPlane(plane);
        double planeOrigin[3];
        plane->GetOrigin(planeOrigin);

        double sliceHeight = planeOrigin[2];
        double current_slice = (sliceHeight - dataOrigin[2]) / dataSpacing;

        horizontalWindSlicer->SetHeight((int) current_slice);
        horizontalWindSlicer->Update();

        windMask->SetInputData(horizontalWindSlicer->GetOutput());

        windMask->Update();
        windGlyphs->Update();
        windGlyphSource->Update();
        windMapper->Update();
        renderer->RemoveActor(horizontalWindVectorActor);

        horizontalWindVectorActor = vtkNew<vtkActor>();
        horizontalWindVectorActor->SetMapper(windMapper);
        horizontalWindVectorActor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
        double position[3];
        horizontalWindVectorActor->GetPosition(position);
        position[2] = sliceHeight;
        horizontalWindVectorActor->SetPosition(position);

        if(windMode == 2){
            horizontalWindVectorActor->SetVisibility(false);
        }

        renderer->AddActor(horizontalWindVectorActor);
        renderer->Render();

        renderWindow->Render();
    }

    void PlaneWidgetUpdated() {
        UpdatePressureSlicer();
        UpdateWindSlice();
        UpdateDivergenceSlicer();

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

    void ToggleCloudRendering() {
        cloudVisible = !cloudVisible;

        cloudActor->SetVisibility(cloudVisible);
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

    void ToggleWindMode() {
        windMode = (windMode + 1) % 3;
        if (windMode == 2) {
            horizontalWindVectorActor->SetVisibility(false);
        } else {
            horizontalWindVectorActor->SetVisibility(true);
        }

    }

    void UpdateDivergenceSlicer() {
        if (!divergenceVisible) return;

        vtkNew<vtkPlane> plane;
        planeWidget->GetPlane(plane);

        horizontalWindDivergenceSliceMapper->SetSlicePlane(plane);

    }

    void ToggleDivergence() {
        divergenceVisible = !divergenceVisible;
        if (!divergenceVisible && pressureVisible) {
            TogglePressure();
        }
        UpdateDivergenceSlicer();

        horizontalWindDivergenceActor->SetVisibility(divergenceVisible);

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
        if (key == "3") {
            ToggleCloudRendering();
        }
        if (key == "4") {
            ToggleWindMode();
        }
        if (key == "5") {
            ToggleDivergence();
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
