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
                    0, 0.7, 0, // green
                    0, 0.6, 0.8,
                    0, 0.4000, 0.6745,
                    0, 0.5000, 0.7745,
                    0, 0.5000, 0.8,
                    0, 0.6000, 0.85,
                    0, 0.6000, 0.9,
                    0, 0.7000, 0.95,
                    0, 0.7000, 1  // Blue
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
    // ----------------------------------------------------------------
    // color initialization based on height
    // ----------------------------------------------------------------
    int elementsx = 1429;//how many x values on grid
    int elementsy = 1556;//how many y
    vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
    colors->SetNumberOfValues(elementsx * elementsy);
    // ----------------------------------------------------------------
    // read in text file with height data to decide on color
    // ----------------------------------------------------------------
    std::cout << getDataPath("/data/height_values.txt") << std::endl;
    std::ifstream myfile(getDataPath("/data/height_values.txt"));
    if (!myfile.is_open()) {
        std::cout << "error opening file" << std::endl;
        return -1;
    }
    for (unsigned int ix = 0; ix < elementsx; ix++) {
        for (unsigned int iy = 0; iy < elementsy; iy++) {
            std::string str;
            std::getline(myfile, str);
            double h = std::stod(str) * 0.0001;
            // Height value color-coded
            colors->SetValue(ix * elementsy + iy, h);
        }
    }
    // ----------------------------------------------------------------
    // Create a lookup table to share between the mapper and the scalar bar
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetScaleToLinear();
    lookupTable->SetNumberOfTableValues(numColors);
    double r, g, b;
    for (int i = 0; i < numColors; i++) {
        double val = min + ((double) i / numColors) * (max - min);
        getColorCorrespondingTovalue(val, r, g, b);
        lookupTable->SetTableValue(i, r, g, b);
    }
    lookupTable->Build();
    // ----------------------------------------------------------------
    // Create a scalar bar actor for the colormap
    // ----------------------------------------------------------------
    vtkSmartPointer<vtkScalarBarActor> legend = vtkSmartPointer<vtkScalarBarActor>::New();
    legend->SetLookupTable(lookupTable);
    legend->SetNumberOfLabels(3);
    legend->SetTitle("height");
    legend->SetVerticalTitleSeparation(6);
    legend->GetPositionCoordinate()->SetValue(0.88, 0.1);
    legend->SetWidth(0.1);
    // ----------------------------------------------------------------
    // Create a triangulated mapper for mesh read from file, color it and link it to the lookup table
    // ----------------------------------------------------------------
    vtkNew<vtkOBJReader> reader;
    reader->SetFileName(getDataPath("/data/heightfield.obj").c_str());
    reader->Update();//read .obj file
    vtkSmartPointer<vtkPolyDataMapper> triangulatedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    triangulatedMapper->SetInputConnection(reader->GetOutputPort());//take mesh as input
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