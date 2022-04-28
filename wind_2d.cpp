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


int main(int, char* []) {

    // ----------------------------------------------------------------
    // read the field data from PolyData.
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName("C:/Users/david/Desktop/SciVis/Clouds/2d_velocity/poly_data.vtp");
    reader->Update();
    vtkSmartPointer<vtkPolyData> v = reader->GetOutput();
   
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
    //densityFilter->SetInputConnection(reader->GetOutputPort());
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

    // ----------------------------------------------------------------
    // Rendering
    // ----------------------------------------------------------------
    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(255,255,255);
    renderer->AddActor(vectors_actor);
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("Top view of Vector Field");

    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindow->Render();
    renderWindowInteractor->Initialize();

    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}