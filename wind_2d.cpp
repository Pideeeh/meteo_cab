#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkGlyph2D.h>
#include <vtkGlyphSource2D.h>
#include <vtkMaskPoints.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkFloatArray.h>
#include <vtkStreamTracer.h>
#include <vtkPolyDataPointSampler.h> 
#include <vtkEvenlySpacedStreamlines2D.h>
#include <vtkPlaneSource.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkButtonWidget.h>
#include <vtkCommand.h>
#include <vtkTexturedButtonRepresentation2D.h>

#include "dataSetDefinitions.h"

int main(int, char* []) {

    vtkNew<vtkNamedColors> colors;

    // ----------------------------------------------------------------
    // read the field data from ImageData.
    // ----------------------------------------------------------------

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(getDataPath("/data/2d_velocity.vti").c_str());
    reader->Update();
    vtkSmartPointer<vtkImageData> v = reader->GetOutput();
   
    // ----------------------------------------------------------------
    // Change spacing (Peter's code), it does not work with polydata
    // ----------------------------------------------------------------
    //double* bounds = v->GetSpacing();
    //bounds[2] *= 0.0001;    // scale the z spacing down
    //v->SetSpacing(bounds);

    v->GetPointData()->SetActiveVectors("2d_velocity");

    // ----------------------------------------------------------------
    // Mask data for nicer visualization
    // ----------------------------------------------------------------
    auto mask = vtkSmartPointer<vtkMaskPoints>::New();
    mask->SetInputData(v);
    mask->RandomModeOn();
    mask->SetRandomModeType(2);
    mask->SetMaximumNumberOfPoints(1000);
    mask->SetOnRatio(50);
    mask->Update();

    // ----------------------------------------------------------------
    // Prepare data for Glyph visualization
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkPolyData> mask_output = mask->GetOutput();
    mask_output->GetPointData()->SetActiveScalars("2d_velocity");

    // ----------------------------------------------------------------
    // Apply Glyphs
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkGlyphSource2D> glyph_source = vtkSmartPointer<vtkGlyphSource2D>::New();
    glyph_source->SetGlyphTypeToArrow();
    glyph_source->FilledOff();

    vtkSmartPointer<vtkGlyph2D> glyph = vtkSmartPointer<vtkGlyph2D>::New();
    glyph->SetInputData(mask_output);
    glyph->SetSourceConnection(glyph_source->GetOutputPort());
    glyph->OrientOn();
    glyph->SetScaleModeToScaleByVector();
    glyph->SetScaleFactor(0.02);
    glyph->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph->GetOutputPort());
    mapper->Update();
    auto vectors_actor = vtkSmartPointer<vtkActor>::New();
    vectors_actor->SetMapper(mapper);
    vectors_actor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());

    
    // ----------------------------------------------------------------
    // Compute Seeds - uniformly distributed in the domain
    // ----------------------------------------------------------------
    v->GetPointData()->SetActiveScalars("2d_velocity");

    double N = 15;
    double bounds[6];
    v->GetBounds(bounds);

    double range[3];
    for (int i = 0; i < 3; ++i) range[i] = bounds[2 * i + 1] - bounds[2 * i];

    vtkNew<vtkPoints> array_seeds;
    for (int y = 0; y <= N; y++) {
        for (int x = 0; x <= N; x++) {
            int idx = y * N + x;
            array_seeds->InsertNextPoint(
                bounds[0] + (x / N) * (range[0]),
                bounds[2] + (y / N) * (range[1]),
                0);
        }
    }

    vtkSmartPointer<vtkPolyData> seeds = vtkSmartPointer<vtkPolyData>::New();
    seeds->SetPoints(array_seeds);

    // ----------------------------------------------------------------
    // Setup streamlines and integrator
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkStreamTracer> tracer = vtkSmartPointer<vtkStreamTracer>::New();
    tracer->SetInputData(v);
    tracer->SetSourceData(seeds);
    tracer->SetMaximumPropagation(1);
    tracer->SetIntegratorTypeToRungeKutta45();
    tracer->SetInitialIntegrationStep(0.01);
    tracer->SetIntegrationDirectionToBoth();
    tracer->Update();

    vtkNew<vtkPolyDataMapper> lines_mapper;
    lines_mapper->SetInputConnection(tracer->GetOutputPort());
    vtkNew<vtkActor> streamlines_actor;
    streamlines_actor->SetMapper(lines_mapper);
    streamlines_actor->VisibilityOn();

    // ----------------------------------------------------------------
    // Rendering
    // ----------------------------------------------------------------
    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(colors->GetColor3d("White").GetData());
    //renderer->AddActor(vectors_actor);
    renderer->AddActor(streamlines_actor);
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("Top view of Vector Field");

    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);
    vtkNew<vtkInteractorStyleTrackballCamera> trackball;
    renderWindowInteractor->SetInteractorStyle(trackball);

    renderWindow->Render();
    renderWindowInteractor->Initialize();

    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}