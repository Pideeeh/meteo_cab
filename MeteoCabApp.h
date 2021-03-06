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
#include "vtkPNGReader.h"
#include "vtkXMLImageDataReader.h"
#include "vtkImageReslice.h"
#include "vtkDoubleArray.h"
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
#include "vtkImageSliceMapper.h"
#include "vtkProperty.h"
#include "vtkOBJReader.h"
#include "vtkPlaneWidget.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMetaImageReader.h"
#include "vtkVolumeProperty.h"
#include "vtkEasyTransfer.hpp"
#include "vtkTransformPolyDataFilter.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkMaskPoints.h"
#include "vtkGlyphSource2D.h"
#include "vtkGradientFilter.h"
#include "vtkInteractorStyleFlight.h"
#include "vtkGlyph2D.h"

#include "utils.h"
#include "dataSetDefinitions.h"
#include "PlaneWidgetInteraction.h"
#include "BetterVtkSlicer.h"
#include "streamline_functions.h"
#include "vtkLegendBoxActor.h"
#include "vtkCubeSource.h"
#include "vtkTextProperty.h"


#define interactorBaseClass vtkInteractorStyleTrackballCamera
#define interactorBaseClass2 vtkInteractorStyleFlight

using namespace std;


class MeteoCabApp : public interactorBaseClass {
    static MeteoCabApp *New() {
        return new MeteoCabApp();
    }

vtkTypeMacro(MeteoCabApp, interactorBaseClass);
    vtkNew<vtkNamedColors> colors;

    /*
     * Pressure Variables
     * */
    bool pressureVisible = false;
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

    vtkSmartPointer<vtkVolume> divergenceActor;
    vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> divergenceRaycastMapper;
    vtkSmartPointer<vtkScalarBarActor> divLegend;

    vtkSmartPointer<vtkActor> horizontalWindVectorActor;
    vtkSmartPointer<vtkImageData> horizontalWindData;
    vtkSmartPointer<vtkImageData> streamlinesData;
    vtkSmartPointer<BetterVtkSlicer> horizontalWindSlicer;
    vtkSmartPointer<vtkPolyDataMapper> windMapper;
    vtkSmartPointer<vtkGlyph2D> windGlyphs;
    vtkSmartPointer<vtkGlyphSource2D> windGlyphSource;
    vtkSmartPointer<vtkMaskPoints> windMask;

    vtkSmartPointer<vtkPolyDataMapper> streamlinesMapper;
    vtkSmartPointer<vtkActor> streamlinesActor;
/*
 * Rendering Variables
 * */
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

    /*
     * Cloud variables
     * */
    bool cliVisible = false;
    vtkSmartPointer<vtkImageData> cliData;
    vtkSmartPointer<vtkVolume> cliActor;
    vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> cliRaycastMapper;

    bool clwVisible = false;
    vtkSmartPointer<vtkImageData> clwData;
    vtkSmartPointer<vtkVolume> clwActor;
    vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> clwRaycastMapper;

    vtkSmartPointer<vtkLegendBoxActor> cli_clwLegend;

    /*
    * Rain variables
    * */
    bool qrVisible = false;
    vtkSmartPointer<vtkImageData> qrData;
    vtkSmartPointer<vtkVolume> qrActor;
    vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> qrRaycastMapper;

    vtkSmartPointer<vtkLegendBoxActor> qrLegend;

    /*
     * LIC
     *
     */

    bool licVisible = false;
    vtkSmartPointer<vtkImageData> licImage;
    vtkSmartPointer<vtkImageMapper3D> licImageMapper;
    vtkSmartPointer<vtkImageActor> licActor;
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
        renderer->AddActor(cliActor);
        renderer->AddActor(clwActor);
        renderer->AddActor(cli_clwLegend);
        renderer->AddActor(qrActor);
        renderer->AddActor(qrLegend);
        renderer->AddActor(divergenceActor);
        renderer->AddActor2D(heightmapLegend);
        renderer->AddActor2D(divLegend);
        renderer->AddActor(licActor);
        renderer->AddActor(horizontalWindVectorActor);
        renderer->AddActor(streamlinesActor);
        renderer->SetBackground(0.7, 0.7, 0.7);
        interactor->SetInteractorStyle(this);
        renderWindow->SetSize(500, 500);
        renderWindow->Render();
        interactor->Start();
    }

    void InitializePressureSlicer() {
        vtkNew<vtkXMLImageDataReader> reader;
        reader->SetFileName(getDataPath("/data/pres_regularly_resampled_downsampled.vti").c_str());
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
        double height_min = 0;//initialization on values of our height data
        double height_max = 10 * 0.244379;
        double numColors = 20;

        // ----------------------------------------------------------------
        // Create a triangulated mapper for mesh read from file, color it and link it to the lookup table
        // ----------------------------------------------------------------
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName(getDataPath("/data/strong_boundaries.obj").c_str());
        reader->Update();//read .obj file
        heightmapDataMapper = vtkNew<vtkPolyDataMapper>();
        heightmapDataMapper->SetInputConnection(reader->GetOutputPort());

        heightmapColors = vtkNew<vtkFloatArray>();
        heightmapColors->SetNumberOfValues(577290);
        for (int i = 0; i < heightmapColors->GetNumberOfValues(); i++) {//go over vertices
            heightmapColors->SetValue(i, 10 * heightmapDataMapper->GetInput()->GetPoint(i)[2]);
            if (heightmapDataMapper->GetInput()->GetPoint(i)[2] == 0) {//sea
                heightmapColors->SetValue(i, NAN);
            }
        }
        // ----------------------------------------------------------------
        // Create a lookup table to share between the mapper and the scalar bar
        // ----------------------------------------------------------------
        heightmapLookupTable = vtkNew<vtkLookupTable>();
        heightmapLookupTable->SetScaleToLinear();
        heightmapLookupTable->SetNumberOfTableValues(numColors);
        heightmapLookupTable->SetNanColor(130. / 255, 167. / 255, 196. / 255, 1);//sea is blue but discrete cutoff
        double r, g, b;
        for (int i = 0; i < numColors; i++) {
            double val = height_min + ((double) i / numColors) * (height_max - height_min);
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
        heightmapLegend->SetTitle("height [km]");
        heightmapLegend->SetVerticalTitleSeparation(6);
        heightmapLegend->GetPositionCoordinate()->SetValue(0.85, 0.1);
        heightmapLegend->SetWidth(0.05);

        heightmapDataMapper->GetInput()->GetPointData()->SetScalars(
                heightmapColors);//color it based on our computations
        heightmapDataMapper->SetLookupTable(heightmapLookupTable);//link lookup table
        heightmapDataMapper->SetScalarRange((height_min), (height_max));
        // ----------------------------------------------------------------
        // Create an actor
        // ----------------------------------------------------------------
        heightmapActor = vtkNew<vtkActor>();
        heightmapActor->SetMapper(heightmapDataMapper);
        heightmapActor->GetProperty()->SetEdgeVisibility(false);
        double origin[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = pressureData->GetOrigin()[2] - 106.81226 * 0.0001;
        heightmapActor->SetPosition(origin);
    }

    void InitializePlaneWidget() {
        planeWidget = vtkNew<vtkPlaneWidget>();
        planeWidget->SetInteractor(interactor);
        planeWidget->SetRepresentationToOutline();

        planeWidgetInteraction = vtkNew<PlaneWidgetInteraction>();
        planeWidgetInteraction->dataSpace = pressureData;
        planeWidgetInteraction->widget = planeWidget;
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

    void InitializeHorizontalWind() {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/2d_wind_full_regularly_resampled_downsampled.vti").c_str());
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
        windGlyphs->SetScaleFactor(0.02);
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


        streamlinesData = vtkNew<vtkImageData>();
        streamlinesData->DeepCopy(horizontalWindData);
        UpdateStreamlines(0);

    }

    void UpdateStreamlines(int streamSlice) {
        if (streamSlice == 0)
            streamSlice = 1;
        streamlinesData->GetPointData()->SetActiveScalars("2d_velocity");

        double spacing[3];
        double origin[3];
        streamlinesData->GetSpacing(spacing);
        streamlinesData->GetOrigin(origin);

        double sliceHeight = streamSlice * spacing[2] + origin[2];

        double N = 15;
        double bounds[6];
        streamlinesData->GetBounds(bounds);
        double range[3];
        for (int i = 0; i < 3; ++i) range[i] = bounds[2 * i + 1] - bounds[2 * i];

        vtkNew<vtkPoints> array_seeds;
        for (int y = 0; y <= N; y++) {
            for (int x = 0; x <= N; x++) {
                array_seeds->InsertNextPoint(
                        bounds[0] + (x / N) * (range[0]),
                        bounds[2] + (y / N) * (range[1]),
                        sliceHeight);
            }
        }

        vtkSmartPointer<vtkPolyData> seeds = vtkSmartPointer<vtkPolyData>::New();
        seeds->SetPoints(array_seeds);


        vtkDataArray *vectors = streamlinesData->GetPointData()->GetVectors();
        streamlinesData->GetPointData()->SetActiveScalars("2d_velocity");
        int dim[3];
        streamlinesData->GetDimensions(dim);
        int z = streamSlice;
        for (int x = 0; x < dim[0]; x++) {
            for (int y = 0; y < dim[1]; y++) {
                auto pixel = static_cast<float *>(streamlinesData->GetScalarPointer(x, y, 0));
                auto kd = vectors->GetTuple(x + y * dim[0] + z * dim[0] * dim[1]);
                pixel[0] = kd[0];
                pixel[1] = kd[1];
                pixel[2] = 0.0;
            }
        }



        // ----------------------------------------------------------------
        // TODO: test code for new streamlines visualization
        vtkNew<vtkPolyData> test1;
        vtkNew<vtkPoints> all_points;
        vtkNew<vtkCellArray> all_lines;
        vtkNew<vtkUnsignedCharArray> all_colors;
        all_colors->SetNumberOfComponents(3);

        for (int i = 0; i < seeds->GetNumberOfPoints(); i++) {
            compute_one_streamline(streamlinesData, seeds->GetPoint(i), test1);
            from_points_to_line(test1->GetPoints(), all_points, all_lines, all_colors);
        }
        test1->SetPoints(all_points);
        test1->SetLines(all_lines);
        test1->GetCellData()->SetScalars(all_colors);

        // ----------------------------------------------------------------
        streamlinesMapper = vtkNew<vtkPolyDataMapper>();
        streamlinesMapper->SetInputData(test1);

        streamlinesActor = vtkNew<vtkActor>();
        streamlinesActor->SetMapper(streamlinesMapper);
        streamlinesActor->VisibilityOn();
        streamlinesActor->GetProperty()->SetLineWidth(1.2);

        streamlinesActor->SetPosition(0, 0, sliceHeight - test1->GetBounds()[4]);

        streamlinesActor->SetVisibility(windMode == 2);

    }

    void ComputeWindDivergence() {

        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/2d_wind_full_regularly_resampled_downsampled.vti").c_str());
        reader->Update();

        vtkNew<vtkGradientFilter> gradientFilter;
        gradientFilter->SetInputConnection(0, reader->GetOutputPort());
        gradientFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "2d_velocity");
        gradientFilter->ComputeDivergenceOn();
        gradientFilter->Update();

        auto div = gradientFilter->GetOutput();
        vtkImageData *input = vtkImageData::SafeDownCast(div);
        input->GetPointData()->SetActiveScalars("Divergence");
        horizontalWindDivergence = input;

        vtkSmartPointer<vtkImageReslice> slicer = vtkNew<vtkImageReslice>();
        slicer->SetInputData(horizontalWindDivergence);

        double origin[3];
        horizontalWindDivergence->GetOrigin(origin);

        double spacing[3];
        horizontalWindDivergence->GetSpacing(spacing);

        int extent[6];
        horizontalWindDivergence->GetExtent(extent);
        origin[0] += spacing[0] * 20;
        origin[1] += spacing[1] * 20;
        extent[1] -= 20;
        extent[3] -= 20;
        extent[5] -= 100;
        slicer->SetOutputExtent(extent);
        slicer->SetOutputOrigin(origin);
        slicer->Update();

        horizontalWindDivergence = slicer->GetOutput();


        double range[2];
        string name = horizontalWindDivergence->GetPointData()->GetArray(0)->GetName();
        horizontalWindDivergence->GetPointData()->SetActiveScalars(name.c_str());
        horizontalWindDivergence->GetPointData()->GetScalars(name.c_str())->GetRange(range);

        divergenceRaycastMapper = vtkNew<vtkOpenGLGPUVolumeRayCastMapper>();
        divergenceRaycastMapper->SetInputData(horizontalWindDivergence);

        auto opacityFunction = GetDivergenceOpacityFunction(range[0], range[1], 0.94, -65, 65);
        auto colorFunction = GetDivergenceColorFunction(range[0], range[1]);

        // assign transfer function to volume properties
        vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        volumeProperty->SetColor(colorFunction);
        volumeProperty->SetScalarOpacity(opacityFunction);

        // ----------------------------------------------------------------
        // Create a lookup table to share between the mapper and the scalar bar
        // ----------------------------------------------------------------
        vtkNew<vtkLookupTable> divLookupTable = vtkNew<vtkLookupTable>();
        divLookupTable->SetScaleToLinear();
        divLookupTable->SetNumberOfTableValues(20);
        for (int i = 0; i < divLookupTable->GetNumberOfTableValues(); i++) {
            double val = range[0] + ((double) i / divLookupTable->GetNumberOfTableValues()) * (range[1] - range[0]);
            double *col = colorFunction->GetColor(val);
            double op = opacityFunction->GetValue(val);
            divLookupTable->SetTableValue(i, col[0], col[1], col[2], op);
        }
        divLookupTable->SetRange(range[0], range[1]);
        divLookupTable->Build();
        divLegend = vtkNew<vtkScalarBarActor>();
        divLegend->SetLookupTable(divLookupTable);
        divLegend->SetNumberOfLabels(3);
        divLegend->SetTitle("divergence");
        divLegend->SetVerticalTitleSeparation(6);
        divLegend->GetPositionCoordinate()->SetValue(0.77, 0.1);
        divLegend->SetWidth(0.05);

        //divergenceRaycastMapper->SetLookupTable(divLookupTable);//link lookup table
        //divergenceRaycastMapper->SetScalarRange(range[0], range[1]);

        // create volume actor and assign mapper and properties
        divergenceActor = vtkNew<vtkVolume>();
        divergenceActor->SetMapper(divergenceRaycastMapper);
        divergenceActor->SetProperty(volumeProperty);
        divergenceActor->SetVisibility(divergenceVisible);
        divLegend->SetVisibility(divergenceVisible);
    }

    void InitializeClouds() {
        //CLI
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/cli_regularly_resampled_downsampled.vti").c_str());
        reader->Update();

        cliData = reader->GetOutput();
        cliData->GetPointData()->SetActiveScalars("cli");

        cliRaycastMapper = vtkNew<vtkOpenGLGPUVolumeRayCastMapper>();
        cliRaycastMapper->SetInputData(cliData);

        double cliRange[2];
        cliData->GetPointData()->GetScalars("cli")->GetRange(cliRange);

        //for cli
        auto cliOpacityFunction = GetSimpleOpacityFunction(cliRange[1], 0.85, cliRange[1] / 25.0);
        auto cliColorFunction = GetSimpleColorFunction(cliRange[1], 1.0, 1.0, 1.0);

        // assign transfer function to volume properties
        vtkSmartPointer<vtkVolumeProperty> cliVolumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        cliVolumeProperty->SetColor(cliColorFunction);
        cliVolumeProperty->SetScalarOpacity(cliOpacityFunction);
        cliVolumeProperty->SetInterpolationTypeToLinear();

        // create volume actor and assign mapper and properties
        cliActor = vtkNew<vtkVolume>();
        cliActor->SetMapper(cliRaycastMapper);
        cliActor->SetProperty(cliVolumeProperty);
        cliActor->SetVisibility(cliVisible);


        //CLW
        reader->SetFileName(getDataPath("/data/clw_regularly_resampled_downsampled.vti").c_str());
        reader->Update();

        clwData = reader->GetOutput();
        clwData->GetPointData()->SetActiveScalars("clw");

        clwRaycastMapper = vtkNew<vtkOpenGLGPUVolumeRayCastMapper>();
        clwRaycastMapper->SetInputData(clwData);

        double clwRange[2];
        clwData->GetPointData()->GetScalars("clw")->GetRange(clwRange);

        //for clw
        auto clwOpacityFunction = GetSimpleOpacityFunction(clwRange[1], 0.9, clwRange[1] / 30.0);
        auto clwColorFunction = GetSimpleColorFunction(clwRange[1], 0.4, 0.4, 0.4);

        // assign transfer function to volume properties
        vtkSmartPointer<vtkVolumeProperty> clwVolumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        clwVolumeProperty->SetColor(clwColorFunction);
        clwVolumeProperty->SetScalarOpacity(clwOpacityFunction);
        clwVolumeProperty->SetInterpolationTypeToLinear();

        // create volume actor and assign mapper and properties
        clwActor = vtkNew<vtkVolume>();
        clwActor->SetMapper(clwRaycastMapper);
        clwActor->SetProperty(clwVolumeProperty);
        clwActor->SetVisibility(clwVisible);


        //add categorical legend
        cli_clwLegend = vtkSmartPointer<vtkLegendBoxActor>::New();

        cli_clwLegend->SetNumberOfEntries(2);

        double color[4];

        vtkNew<vtkCubeSource> legendBox;
        legendBox->Update();
        color[0] = 1.0;
        color[1] = 1.0;
        color[2] = 1.0;
        color[3] = 1.0;
        cli_clwLegend->SetEntry(0, legendBox->GetOutput(), "Cloud ice content", color);

        color[0] = 0.4;
        color[1] = 0.4;
        color[2] = 0.4;
        color[3] = 1.0;
        cli_clwLegend->SetEntry(1, legendBox->GetOutput(), "Cloud water content", color);

        // place legend
        cli_clwLegend->GetPositionCoordinate()->SetValue(0.02, 0.8);
        cli_clwLegend->GetPosition2Coordinate()->SetValue(0.16, 0.12);

        cli_clwLegend->UseBackgroundOn();
        double background[4];
        background[0] = 0.0;
        background[1] = 0.0;
        background[2] = 0.0;
        background[3] = 1.0;
        //colors->GetColor("SlateGray", background);
        cli_clwLegend->SetBackgroundColor(background);
        cli_clwLegend->SetVisibility(qrVisible);
    }

    void InitializeRain() {
        //qr aka rain
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(getDataPath("/data/qr_regularly_resampled_downsampled.vti").c_str());
        reader->Update();

        qrData = reader->GetOutput();
        qrData->GetPointData()->SetActiveScalars("qr");

        qrRaycastMapper = vtkNew<vtkOpenGLGPUVolumeRayCastMapper>();
        qrRaycastMapper->SetInputData(qrData);

        double qrRange[2];
        qrData->GetPointData()->GetScalars("qr")->GetRange(qrRange);

        //for qr
        auto qrOpacityFunction = GetSimpleOpacityFunction(qrRange[1], 0.95, qrRange[1] / 50.0);
        auto qrColorFunction = GetSimpleColorFunction(qrRange[1], 0.4, 0.6, 0.99);

        // assign transfer function to volume properties
        vtkSmartPointer<vtkVolumeProperty> qrVolumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        qrVolumeProperty->SetColor(qrColorFunction);
        qrVolumeProperty->SetScalarOpacity(qrOpacityFunction);
        qrVolumeProperty->SetInterpolationTypeToLinear();

        // create volume actor and assign mapper and properties
        qrActor = vtkNew<vtkVolume>();
        qrActor->SetMapper(qrRaycastMapper);
        qrActor->SetProperty(qrVolumeProperty);
        qrActor->SetVisibility(qrVisible);

        //add categorical legend
        qrLegend = vtkSmartPointer<vtkLegendBoxActor>::New();

        qrLegend->SetNumberOfEntries(1);

        double color[4];

        vtkNew<vtkCubeSource> legendBox;
        legendBox->Update();
        color[0] = 0.4;
        color[1] = 0.6;
        color[2] = 0.99;
        color[3] = 1.0;
        qrLegend->SetEntry(0, legendBox->GetOutput(), "Rain mixing ratio   ", color);

        // place legend
        qrLegend->GetPositionCoordinate()->SetValue(0.02, 0.73);
        qrLegend->GetPosition2Coordinate()->SetValue(0.16, 0.06);

        qrLegend->UseBackgroundOn();
        double background[4];
        background[0] = 0.0;
        background[1] = 0.0;
        background[2] = 0.0;
        background[3] = 1.0;
        qrLegend->SetBackgroundColor(background);
        qrLegend->SetVisibility(qrVisible);
    }


    void Launch() {
        InitializePressureSlicer();
        InitializeHeightmap();

        InitializeHorizontalWind();
        ComputeWindDivergence();
        InitializeClouds();
        InitializeRain();
        InitializeLICRendering();
        StartVisualization();
    }

    void UpdatePressureSlicer() {

        int dimensions[3];
        pressureData->GetDimensions(dimensions);
        double bounds[6];
        pressureData->GetBounds(bounds);

        double dataOrigin[3];
        pressureData->GetOrigin(dataOrigin);

        double dataSpacing = pressureData->GetSpacing()[2];
        double origin[3];
        pressureReslicer->GetOutputOrigin(origin);

        vtkNew<vtkPlane> plane;
        planeWidget->GetPlane(plane);

        double planeOrigin[3];
        plane->GetOrigin(planeOrigin);

        double sliceHeight = std::min(bounds[5],
                                      std::max(planeOrigin[2], bounds[4]));//avoid exception when user goes too far
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
        pressureLegend->SetTitle("pressure [Pa]");
        pressureLegend->SetVerticalTitleSeparation(6);
        pressureLegend->GetPositionCoordinate()->SetValue(0.93, 0.1);
        pressureLegend->SetWidth(0.05);
        pressureLegend->SetVisibility(pressureVisible);

        renderer->AddActor2D(pressureLegend);

        pressureContourFilter->GenerateValues(10, pressure_value_lower, pressure_value_upper);
        pressureContourFilter->Update();
        pressureContourMapper->Update();
        pressureSliceActor->Update();

        double *pos = contourActor->GetPosition();
        contourActor->SetPosition(pos);

    }

    void UpdateWindSlice() {
        if (windMode == 0) return;
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

        double bounds[6];
        horizontalWindData->GetBounds(bounds);
        double sliceHeight = std::min(bounds[5],
                                      std::max(planeOrigin[2], bounds[4]));//avoid exception when user goes too far
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

        if (windMode == 2) {
            horizontalWindVectorActor->SetVisibility(false);
        }

        renderer->RemoveActor(streamlinesActor);
        UpdateStreamlines((int) current_slice);
        renderer->AddActor(streamlinesActor);

        renderer->AddActor(horizontalWindVectorActor);
        renderer->Render();

        renderWindow->Render();
    }

    void PlaneWidgetUpdated() {
        UpdatePressureSlicer();
        UpdateWindSlice();

        UpdateLICSlice(false);

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
        cliVisible = !cliVisible;

        cliActor->SetVisibility(cliVisible);

        clwVisible = !clwVisible;

        clwActor->SetVisibility(clwVisible);

        cli_clwLegend->SetVisibility(cliVisible);

    }

    void ToggleRainRendering() {
        qrVisible = !qrVisible;

        qrActor->SetVisibility(qrVisible);

        qrLegend->SetVisibility(qrVisible);

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
        interactorBaseClass::OnLeftButtonUp();
    }

    void ToggleWindMode() {
        windMode = (windMode + 1) % 3;
        if (windMode == 2) {
            horizontalWindVectorActor->SetVisibility(false);
            streamlinesActor->SetVisibility(true);
        } else if (windMode == 1) {
            horizontalWindVectorActor->SetVisibility(true);
            streamlinesActor->SetVisibility(false);
        } else {
            horizontalWindVectorActor->SetVisibility(false);
            streamlinesActor->SetVisibility(false);
        }

    }

    void ToggleDivergence() {
        divergenceVisible = !divergenceVisible;
        divergenceActor->SetVisibility(divergenceVisible);
        divLegend->SetVisibility(divergenceVisible);
    }

    void ToggleLICRendering() {
        licVisible = !licVisible;
        licActor->SetVisibility(licVisible);

    }

    void InitializeLICRendering() {
        UpdateLICSlice(true);

    }

    void UpdateLICSlice(bool init) {
        licImage = vtkNew<vtkImageData>();
        licImageMapper = vtkNew<vtkImageSliceMapper>();


        int selectedSlice = 0;
        if (!init) {
            vtkNew<vtkPlane> plane;
            planeWidget->GetPlane(plane);
            double planeOrigin[3];
            plane->GetOrigin(planeOrigin);

            double dataOrigin[3];
            pressureData->GetOrigin(dataOrigin);

            double dataSpacing = pressureData->GetSpacing()[2];
            double bounds[6];
            pressureData->GetBounds(bounds);
            double sliceHeight = std::min(bounds[5],
                                          std::max(planeOrigin[2], bounds[4]));//avoid exception when user goes too far
            double currentSlice = (sliceHeight - dataOrigin[2]) / dataSpacing;

            selectedSlice = currentSlice;
        }

        vtkNew<vtkXMLImageDataReader> reader;
        reader->SetFileName(getDataPath("/data/LIC/lic" + std::to_string(selectedSlice) + ".vti").c_str());
        reader->Update();
        licImage = reader->GetOutput();

        int dimensions[3];
        licImage->GetDimensions(dimensions);
        double dataBounds[6];
        double dataSpacing[3];
        pressureData->GetBounds(dataBounds);
        pressureData->GetSpacing(dataSpacing);
        double spacing[3];
        spacing[0] = (dataBounds[1] - dataBounds[0]) / (double) dimensions[0];
        spacing[1] = (dataBounds[3] - dataBounds[2]) / (double) dimensions[1];
        spacing[2] = 0;

        double origin[3];
        pressureData->GetOrigin(origin);
        origin[2] += dataSpacing[2] * selectedSlice;
        licImage->SetSpacing(spacing);
        licImage->SetOrigin(origin);

        licImageMapper->SetInputData(licImage);

        if(!init){
            renderer->RemoveActor(licActor);
        }
        licActor = vtkNew<vtkImageActor>();
        licActor->SetMapper(licImageMapper);
        licActor->SetVisibility(licVisible);
        licActor->Update();
        if(!init){
            renderer->AddActor(licActor);
        }
    }

    void OnKeyPress() override {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        string key = rwi->GetKeySym();

        if (key == "1") {
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

        if (key == "6") {
            ToggleRainRendering();
        }
        if (key == "7") {
            ToggleLICRendering();
        }

        if (key == "h") {
            planeWidgetInteraction->ToggleOrientation();
            PlaneWidgetUpdated();
        }

        if (pressureVisible) {
            HandlePressureInputs(key);
        }

        renderWindow->Render();


        // Forward events
        interactorBaseClass::OnKeyPress();
    }

    double pressure_value_lower = 75000.0;
    double pressure_value_upper = 100000.0;


};


#endif //VIS2021_METEOCABAPP_H
