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
#include "vtkImageCast.h"

// LIC
#include "vtkImageMapper.h"
#include "vtkActor2D.h"
#include <vtkCharArray.h>
#include <vtkImageGaussianSmooth.h>

void create_noisy_img(int h_res, int v_res, vtkSmartPointer<vtkImageData> image) {
    vtkNew<vtkCharArray> colors;
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(h_res * v_res * 1);
    colors->SetName("colors");

    // Best random integer generator according to stack overflow
    std::random_device rd_device;
    std::mt19937 mt19937(rd_device());
    std::uniform_int_distribution<int> distrib(0, 255);

    image->SetDimensions(h_res, v_res, 1);
    image->SetOrigin(.5, .5, 0);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    // Construction of a random image according to vtk tutorial https://kitware.github.io/vtk-examples/site/Cxx/Images/BorderPixelSize/
    for (unsigned int x = 0; x < h_res; x++)
    {
        for (unsigned int y = 0; y < v_res; y++)
        {
            int grey_value = distrib(mt19937);

            unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
            pixel[0] = static_cast<unsigned char>(grey_value);
            pixel[1] = static_cast<unsigned char>(grey_value);
            pixel[2] = static_cast<unsigned char>(grey_value);
        }
    }
    image->Modified();
}
void compute_one_streamline(vtkSmartPointer<vtkImageData> v, double *point, vtkSmartPointer<vtkPolyData> result) {

    // Allocate one seed
    vtkNew<vtkPolyData> seeds_LIC;
    vtkNew<vtkPoints> one_seed;
    one_seed->InsertNextPoint(point);
    seeds_LIC->SetPoints(one_seed);

    // Compute streamline
    vtkNew<vtkStreamTracer> tracer_LIC;
    tracer_LIC->SetInputData(v);
    tracer_LIC->SetSourceData(seeds_LIC);
    tracer_LIC->SetMaximumPropagation(1);
    tracer_LIC->SetIntegratorTypeToRungeKutta45();
    tracer_LIC->SetInitialIntegrationStep(0.01);
    tracer_LIC->SetIntegrationDirectionToBoth();
    tracer_LIC->Update();

    // return streamline points
    vtkSmartPointer<vtkPoints> points = tracer_LIC->GetOutput()->GetPoints();
    result->SetPoints(points);
}

void weighted_color(vtkSmartPointer<vtkPolyData> polyData, vtkSmartPointer<vtkImageData> noisy_image, double *bounds, int *dims, double *range, unsigned char *color) {
        polyData->GetNumberOfPoints();
        vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
        int n_points = points->GetNumberOfPoints(), counter = 0;

        int buffer[3];
        buffer[0] = buffer[1] = buffer[2] = 0;
        unsigned char* pixel;
        for (unsigned int i = 0; i < n_points; i ++)
        {
            double* point = points->GetPoint(i);

            // Convert point to image space
            int x = (point[0] - bounds[0]) * (dims[0] / range[0]) - 0.5;
            int y = (point[1] - bounds[2]) * (dims[1] / range[1]) - 0.5;

            if (x > 0 && x < dims[0] && y>0 && y < dims[1]) {
                // std::cout << x << " " << y << " " << std::endl;
                // Get pixel color
                pixel = static_cast<unsigned char*>(noisy_image->GetScalarPointer(x, y, 0));
                buffer[0] += pixel[0];
                buffer[1] += pixel[1];
                buffer[2] += pixel[2];
                counter++;
            }
        }

        if(counter > 0){
            //std::cout << "color: " << buffer[0] << " " << buffer[1] << " " << buffer[2] << std::endl;
            //std::cout << "counter: " << counter << std::endl;
            color[0] = buffer[0] / counter;
            color[1] = buffer[1] / counter;
            color[2] = buffer[2] / counter;
        }
}

void normalize_vecotrs(vtkSmartPointer<vtkImageData> data) {
    int dims[3];
    data->GetDimensions(dims);

    for (auto z = 0; z < dims[2]; z++) {
        for (auto y = 0; y < dims[1]; y++) {
            for (auto x = 0; x < dims[0]; x++) {
                int idx = z * dims[0] * dims[1] + y * dims[0] + x;
                auto vector = static_cast<float*>(data->GetScalarPointer(x, y, z));
                double norm = sqrt(vector[0] * vector[0] + vector[1] * vector[1]); // 2d vector
                vector[0] = vector[0] / norm;
                vector[1] = vector[1] / norm;
                vector[2] = 0;
            }
        }
    }
}
void set_texture_color(vtkSmartPointer<vtkImageData> image, vtkSmartPointer<vtkPolyData> polyData, double* bounds, int* dims, double* range, unsigned char *color) {
    polyData->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    int n_points = points->GetNumberOfPoints();

    for (int i = 0; i < n_points; i++) {
        double* point = points->GetPoint(i);

        int x = (point[0] - bounds[0]) * (dims[0] / range[0]) - 0.5;
        int y = (point[1] - bounds[2]) * (dims[1] / range[1]) - 0.5;
        if (x > 0 && x < dims[0] && y>0 && y < dims[1]) {
            unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
            pixel[0] = static_cast<unsigned char>(color[0]);
            pixel[1] = static_cast<unsigned char>(color[1]);
            pixel[2] = static_cast<unsigned char>(color[2]);
        }
    }
    image->Modified();
}

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
    // ----------------------------------------------------------------
    // Visualize glyphs
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
   
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
    // ----------------------------------------------------------------
    // Visualize streamlines
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    
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
    tracer->SetMaximumPropagation(0.5);
    tracer->SetIntegratorTypeToRungeKutta45();
    tracer->SetInitialIntegrationStep(0.0000000001);
    tracer->SetIntegrationDirectionToBoth();
    tracer->Update();

    vtkNew<vtkPolyDataMapper> lines_mapper;
    lines_mapper->SetInputConnection(tracer->GetOutputPort());
    vtkNew<vtkActor> streamlines_actor;
    streamlines_actor->SetMapper(lines_mapper);
    streamlines_actor->VisibilityOn();

    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    // Line integral convolution
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    
    // ----------------------------------------------------------------
    // Create noise image
    // ----------------------------------------------------------------
    // Construction of a random image according to vtk tutorial https://kitware.github.io/vtk-examples/site/Cxx/Images/BorderPixelSize/

    int dims[3];
    v->GetDimensions(dims);
    dims[0] = dims[1] = 256;
    vtkNew<vtkImageData> noisy_image;
    create_noisy_img(dims[0], dims[1], noisy_image);

    // ----------------------------------------------------------------
    // Map image pixels to streamlines seeds
    // ----------------------------------------------------------------
    vtkNew<vtkPoints> array_seeds_LIC;
    for (int y = 0; y <= dims[0]; y++) {
        for (int x = 0; x <= dims[1]; x++) {
            int idx = y * dims[0] + x;
            array_seeds_LIC->InsertNextPoint(
                bounds[0] + ((x + 0.5) / dims[0]) * (range[0]), // ranges computed above for the streamlines
                bounds[2] + ((y + 0.5) / dims[1]) * (range[1]),
                bounds[4]);
        }
    }


    // ----------------------------------------------------------------
    // Construction of LIC texture from vtkStreamTracer. I compute one streamline at a time and find its color.
    // ----------------------------------------------------------------
    
    // Create result image
    vtkNew<vtkImageData> LIC;
    LIC->SetDimensions(dims[0], dims[1], 1);
    LIC->SetOrigin(.5, .5, 0);
    LIC->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    // Create image colors
    vtkNew<vtkCharArray> colors_LIC;
    colors_LIC->SetNumberOfComponents(3);
    colors_LIC->SetNumberOfTuples(dims[0] * dims[1] * 1);
    colors_LIC->SetName("colors");

    int iterations = array_seeds_LIC->GetNumberOfPoints();
    vtkNew<vtkPolyData> streamlines_LIC;

    //v->GetPointData()->SetActiveScalars("2d_velocity");
    //normalize_vecotrs(v);
    //v->GetPointData()->SetActiveVectors("2d_velocity");
    for(int i=0; i<iterations; i++){

        // Compute streamline
        double *point = array_seeds_LIC->GetPoint(i);
        compute_one_streamline(v, point, streamlines_LIC);

        // Compute weighted color of the line
        unsigned char color[3];
        weighted_color(streamlines_LIC, noisy_image, bounds, dims, range, color);
        //std::cout << "color: " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << std::endl;

        // Color LIC texture
        set_texture_color(LIC, streamlines_LIC, bounds, dims, range, color);
        
        std::cout << "iteration: " << i << " of " << iterations << std::endl;
    }

    // ----------------------------------------------------------------
    // Visualize results
    // ----------------------------------------------------------------
    vtkNew<vtkImageMapper> img_map;
    img_map->SetInputData(LIC);
    img_map->SetColorWindow(255); // width of the color range to map to
    img_map->SetColorLevel(127.5); // center of the color range to map to

    vtkNew<vtkActor2D> image_actor;
    image_actor->SetMapper(img_map);

    // ----------------------------------------------------------------
    // Rendering
    // ----------------------------------------------------------------
    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(colors->GetColor3d("White").GetData());
    //renderer->AddActor(vectors_actor);
    //renderer->AddActor(streamlines_actor);
    renderer->AddActor(image_actor);

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