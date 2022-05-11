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
    static const int numColorNodes = 9;
    double color[numColorNodes][3] =
            {
                    0.1, 0.5, 0.1, // green
                    0.1, 0.5, 0.5,
                    0.1, 0.5, 0.8,
                    0.2, 0.8, 1,  // Blue
                    0.3, 0.9, 1,
                    0.8, 0.95,1,
                    0.98,0.98,1,
                    1,1,1,
                    1,1,1//white
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
    reader->SetFileName(getDataPath("/data/downsampled_terrain.obj").c_str());
    reader->Update();//read .obj file
    vtkSmartPointer<vtkPolyDataMapper> triangulatedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    triangulatedMapper->SetInputConnection(reader->GetOutputPort());//take mesh as input
    //color the mesh
    vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
    colors->SetNumberOfValues(229972);
    for (int i = 0; i < 229972; i++) {//go over vertices
        colors->SetValue(i, triangulatedMapper->GetInput()->GetPoint(i)[2]);
    }
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetScaleToLinear();
    lookupTable->SetNumberOfTableValues(numColors);
    double r, g, b;
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