#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkOBJReader.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include "dataSetDefinitions.h"

double min = 0;//initialization on values of our height data
double max = 0.244379;
static const double numColors = 20;

void getColorCorrespondingTovalue(double val, double &r, double &g, double &b) {
    double range = max - min;
    static const int numColorNodes = 7;
    double color[numColorNodes][3] =
            {
                    103./255,130./255,113./255,//green
                    137./255,166./255,148./255,//a bit lighter green
                    166./255,182./255,155./255,//more brown-ish green
                    185./255, 180./255, 150./255,//brown
                    185./255+10./255, 180./255 + 10. / 255, 150./255 + 10. / 255,//lighter brown
                    //197./255,200./255,176./255,//light brown
                    1,1,1,//white
                    1,1,1,//white
            };

    for (int i = 0; i < (numColorNodes - 1); i++) {
        double currFloor = min + ((double) i / (numColorNodes - 1)) * range;
        double currCeil = min + ((double) (i + 1) / (numColorNodes - 1)) * range;

        if ((val >= currFloor) && (val <= currCeil)) {
            double currFraction = (val - currFloor) / (currCeil - currFloor);
            r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
            g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
            b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
        }
    }
}


int main(int, char *[]) {
    // Create a triangulated mapper for mesh read from file, color it based on height
    vtkNew<vtkOBJReader> reader;
    reader->SetFileName(getDataPath("/data/strong_boundaries.obj").c_str());
    reader->Update();//read .obj file
    vtkSmartPointer<vtkPolyDataMapper> triangulatedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    triangulatedMapper->SetInputConnection(reader->GetOutputPort());//take mesh as input
    //color the mesh
    vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
    colors->SetNumberOfValues(577290);
    //colors->SetNumberOfValues(280551);
    //colors->SetNumberOfValues(2223524);
    //colors->SetNumberOfValues(229972);
    for (int i = 0; i < colors->GetNumberOfValues(); i++) {//go over vertices
        colors->SetValue(i, triangulatedMapper->GetInput()->GetPoint(i)[2]);
        if (colors->GetValue(i) == 0) {//sea
            colors->SetValue(i, NAN);
        }
    }
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetScaleToLinear();
    lookupTable->SetNumberOfTableValues(numColors);
    double r, g, b;
    lookupTable->SetNanColor(130. / 255, 167. / 255, 196. / 255, 1);//sea is blue but discrete cutoff
    for (int i = 0; i < numColors; i++) {
        double val = min + ((double) i / numColors) * (max - min);
        getColorCorrespondingTovalue(val, r, g, b);
        lookupTable->SetTableValue(i, r, g, b);
    }
    lookupTable->SetTableRange(min, max);
    lookupTable->Build();
    vtkSmartPointer<vtkScalarBarActor> legend = vtkSmartPointer<vtkScalarBarActor>::New();
    legend->SetLookupTable(lookupTable);
    legend->SetNumberOfLabels(3);
    legend->SetTitle("height");
    legend->SetVerticalTitleSeparation(6);
    legend->GetPositionCoordinate()->SetValue(0.88, 0.1);
    legend->SetWidth(0.1);
    triangulatedMapper->GetInput()->GetPointData()->SetScalars(colors);//color it based on our computations
    triangulatedMapper->SetLookupTable(lookupTable);//link lookup table
    triangulatedMapper->SetScalarRange(min, max);
    // ----------------------------------------------------------------
    // Create an actor
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkActor> triangulatedActor = vtkSmartPointer<vtkActor>::New();
    triangulatedActor->SetMapper(triangulatedMapper);
    triangulatedActor->GetProperty()->SetEdgeVisibility(false);
    // ----------------------------------------------------------------
    // Create a renderer, render window, and interactor
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1, 1, 1);
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> inter = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    inter->SetRenderWindow(renderWindow);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    inter->SetInteractorStyle(style);
    // ----------------------------------------------------------------
    // Add the actors to the scene
    // ----------------------------------------------------------------
    renderer->AddActor(triangulatedActor);
    renderer->AddActor2D(legend);
    // ----------------------------------------------------------------
    // Render and interact
    // ----------------------------------------------------------------
    renderWindow->Render();
    inter->Start();
    return EXIT_SUCCESS;
}