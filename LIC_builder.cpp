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
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include "dataSetDefinitions.h"
#include "vtkImageCast.h"
#include "vtkImageMapper.h"
#include "vtkActor2D.h"
#include <vtkCharArray.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPNGWriter.h>

void interpolate_velocity(vtkSmartPointer<vtkImageData>& data, double x_location, double y_location, double& u, double& v) {
    double bounds[6];
    data->GetBounds(bounds);
    int dims[3];
    data->GetDimensions(dims);

    x_location = std::min(bounds[1], std::max(x_location, bounds[0]));
    y_location = std::min(bounds[3], std::max(y_location, bounds[2]));

    double between_index_x = (x_location - bounds[0]) / (bounds[1] - bounds[0]) * (dims[0] - 1);
    double between_index_y = (y_location - bounds[2]) / (bounds[3] - bounds[2]) * (dims[1] - 1);

    int lower_index_x = (int)std::floor(between_index_x);
    int upper_index_x = (int)std::ceil(between_index_x);
    int lower_index_y = (int)std::floor(between_index_y);
    int upper_index_y = (int)std::ceil(between_index_y);

    float* vector_low_low = static_cast<float*>(data->GetScalarPointer(lower_index_x, lower_index_y, 0));
    float* vector_low_upper = static_cast<float*>(data->GetScalarPointer(lower_index_x, upper_index_y, 0));
    float* vector_upper_low = static_cast<float*>(data->GetScalarPointer(upper_index_x, lower_index_y, 0));
    float* vector_upper_upper = static_cast<float*>(data->GetScalarPointer(upper_index_x, upper_index_y, 0));

    double t_x = between_index_x - (double)lower_index_x; //influence of upper for x
    double t_y = between_index_y - (double)lower_index_y; //influence of upper for x

    u = t_x * t_y * vector_upper_upper[0] + t_x * (1.0 - t_y) * vector_upper_low[0] + (1.0 - t_x) * t_y * vector_low_upper[0] + (1.0 - t_x) * (1.0 - t_y) * vector_low_low[0];
    v = t_x * t_y * vector_upper_upper[1] + t_x * (1.0 - t_y) * vector_upper_low[1] + (1.0 - t_x) * t_y * vector_low_upper[1] + (1.0 - t_x) * (1.0 - t_y) * vector_low_low[1];
}


void runge_kutta_step(vtkSmartPointer<vtkImageData>& data, double x_i, double y_i, double s, double& x_i1, double& y_i1) {
    double bounds[6];
    data->GetBounds(bounds);

    if (std::min(bounds[1], std::max(x_i1, bounds[0])) != x_i) {
        x_i1 = x_i;
        y_i1 = y_i;
        return;
    }
    if (std::min(bounds[3], std::max(y_i1, bounds[2])) != y_i) {
        x_i1 = x_i;
        y_i1 = y_i;
        return;
    }

    double x_v1, y_v1, x_v2, y_v2, x_v3, y_v3, x_v4, y_v4;

    interpolate_velocity(data, x_i, y_i, x_v1, y_v1);
    double length = std::sqrt(x_v1 * x_v1 + y_v1 * y_v1);
    x_v1 /= length;
    y_v1 /= length;
    interpolate_velocity(data, x_i + s / 2.0 * x_v1, y_i + s / 2.0 * y_v1, x_v2, y_v2);
    length = std::sqrt(x_v2 * x_v2 + y_v2 * y_v2);
    x_v2 /= length;
    y_v2 /= length;
    interpolate_velocity(data, x_i + s / 2.0 * x_v2, y_i + s / 2.0 * y_v2, x_v3, y_v3);
    length = std::sqrt(x_v3 * x_v3 + y_v3 * y_v3);
    x_v3 /= length;
    y_v3 /= length;
    interpolate_velocity(data, x_i + s * x_v3, y_i + s * y_v3, x_v4, y_v4);
    length = std::sqrt(x_v4 * x_v4 + y_v4 * y_v4);
    x_v4 /= length;
    y_v4 /= length;
    double step_x = (x_v1 / 6.0 + x_v2 / 3.0 + x_v3 / 3.0 + x_v4 / 6.0);
    double step_y = (y_v1 / 6.0 + y_v2 / 3.0 + y_v3 / 3.0 + y_v4 / 6.0);
    length = std::sqrt(step_x * step_x + step_y * step_y);
    //length = 1.0;
    x_i1 = x_i + s * step_x / length;
    y_i1 = y_i + s * step_y / length;
}


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

void compute_one_streamline(vtkSmartPointer<vtkImageData> v, double* point, int n_steps, double s, vtkSmartPointer<vtkPolyData> result) {

    vtkNew<vtkPoints> line_points;

    point[2] = 0.3;

    double old_point[3];
    old_point[0] = point[0];
    old_point[1] = point[1];
    old_point[2] = 0.3;

    for (int i = 0; i < n_steps; i++) {
        double new_point[3];
        runge_kutta_step(v, old_point[0], old_point[1], s, new_point[0], new_point[1]);
        new_point[2] = 0.3;
        line_points->InsertNextPoint(new_point);

        old_point[0] = new_point[0];
        old_point[1] = new_point[1];
    }
    result->SetPoints(line_points);
}

void weighted_color(vtkSmartPointer<vtkPolyData> polyData, vtkSmartPointer<vtkImageData> noisy_image, double* bounds, int* dims, double* range, unsigned char* color) {
    polyData->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    int n_points = points->GetNumberOfPoints(), counter = 0;

    int buffer[3];
    buffer[0] = buffer[1] = buffer[2] = 0;
    unsigned char* pixel;
    for (unsigned int i = 0; i < n_points; i++)
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

    if (counter > 0) {
        //std::cout << "color: " << buffer[0] << " " << buffer[1] << " " << buffer[2] << std::endl;
        //std::cout << "counter: " << counter << std::endl;
        color[0] = buffer[0] / counter;
        color[1] = buffer[1] / counter;
        color[2] = buffer[2] / counter;
    }
}

void set_texture_color(vtkSmartPointer<vtkImageData> image, int* point, double* bounds, int* dims, double* range, unsigned char* color) {
    int x = point[0];
    int y = point[1];
    unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
    pixel[0] = static_cast<unsigned char>(color[0]);
    pixel[1] = static_cast<unsigned char>(color[1]);
    pixel[2] = static_cast<unsigned char>(color[2]);
    image->Modified();
}

void from_points_to_line(
    vtkSmartPointer<vtkPoints> line_points,
    vtkSmartPointer<vtkPoints> all_points,
    vtkSmartPointer<vtkCellArray> all_lines,
    vtkSmartPointer< vtkUnsignedCharArray> all_colors
) {
    int n_points = line_points->GetNumberOfPoints();
    int offset = all_points->GetNumberOfPoints();

    vtkNew<vtkNamedColors> colors;

    for (int i = 0; i < n_points - 1; i++) {
        // create a segment between each set of point
        double* point = line_points->GetPoint(i);
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, offset + i);
        line->GetPointIds()->SetId(1, offset + i + 1);

        // update set of all points, colors, and lines
        all_points->InsertNextPoint(point);
        all_lines->InsertNextCell(line);
        all_colors->InsertNextTuple(colors->GetColor3d("Blue").GetData());
    }

    // Add the last point to the array
    double* point = line_points->GetPoint(n_points - 1);
    all_points->InsertNextPoint(point);
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
    // Prepare vectors
    // ----------------------------------------------------------------
    vtkDataArray* vectors = v->GetPointData()->GetVectors();
    v->GetPointData()->SetActiveScalars("2d_velocity");
    int dim[3];
    v->GetDimensions(dim);
    for (int x = 0; x < dim[0]; x++) {
        for (int y = 0; y < dim[1]; y++) {
            auto pixel = static_cast<float*>(v->GetScalarPointer(x, y, 0));
            auto kd = vectors->GetTuple(x + y * dim[0]);
            pixel[0] = kd[0];
            pixel[1] = kd[1];
            pixel[2] = 0.0;
        }
    }

    for (int i = 0; i < 3; i++) {

        double bounds[6];
        v->GetBounds(bounds);
        double range[3];
        for (int i = 0; i < 3; ++i) range[i] = bounds[2 * i + 1] - bounds[2 * i];

        // ----------------------------------------------------------------
        // Create noise image
        // ----------------------------------------------------------------
        // Construction of a random image according to vtk tutorial https://kitware.github.io/vtk-examples/site/Cxx/Images/BorderPixelSize/

        int img_dim[2];
        img_dim[0] = 50;
        img_dim[1] = (img_dim[0] * range[1]) / range[0];
        vtkNew<vtkImageData> noisy_image;
        noisy_image->SetSpacing(v->GetSpacing());
        create_noisy_img(img_dim[0], img_dim[1], noisy_image);

        // ----------------------------------------------------------------
        // Construction of LIC texture from vtkStreamTracer. I compute one streamline at a time and find its color.
        // ----------------------------------------------------------------

        vtkNew<vtkImageData> LIC;
        LIC->SetSpacing(v->GetSpacing());
        LIC->SetDimensions(img_dim[0], img_dim[1], 1);
        LIC->SetOrigin(.5, .5, 0);
        LIC->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

        vtkNew<vtkCharArray> colors_LIC;
        colors_LIC->SetNumberOfComponents(3);
        colors_LIC->SetNumberOfTuples(img_dim[0] * img_dim[1] * 1);
        colors_LIC->SetName("colors");

        vtkNew<vtkPolyData> streamlines_LIC;
        int n_steps = 15;
        double s = 0.01;

        for (int y = 0; y < img_dim[1]; y++) {
            for (int x = 0; x < img_dim[0]; x++) {

                // find the pixel
                int pixel[3];
                pixel[0] = x;
                pixel[1] = y;
                pixel[2] = 0;

                // Compute streamline forward
                double point[3];
                point[0] = bounds[0] + ((x + 0.5) / img_dim[0]) * (range[0]);
                point[1] = bounds[2] + ((y + 0.5) / img_dim[1]) * (range[1]);
                point[2] = bounds[4];

                compute_one_streamline(v, point, n_steps, s, streamlines_LIC);

                // Compute weighted color of the line forward
                unsigned char color_forward[3];
                weighted_color(streamlines_LIC, noisy_image, bounds, img_dim, range, color_forward);

                // Compute streamline backward
                compute_one_streamline(v, point, n_steps, -s, streamlines_LIC);

                // Compute weighted color of the line backward
                unsigned char color_backward[3];
                weighted_color(streamlines_LIC, noisy_image, bounds, img_dim, range, color_backward);

                // Color LIC texture
                color_forward[0] = (color_forward[0] + color_backward[0]) / 2;
                color_forward[1] = (color_forward[1] + color_backward[1]) / 2;
                color_forward[2] = (color_forward[2] + color_backward[2]) / 2;
                set_texture_color(LIC, pixel, bounds, img_dim, range, color_forward);

                // Counter for debugging
                int i = x + y * img_dim[0];
                std::cout << "iteration: " << i << " of " << img_dim[0] * img_dim[1] << std::endl;
            }
        }

        // ----------------------------------------------------------------
        // Writing files
        // ----------------------------------------------------------------
        std::string path_vti = getDataPath("/data/test" + std::to_string(i) + ".vti");
        vtkNew<vtkXMLImageDataWriter> writer_vti;
        writer_vti->SetFileName(path_vti.c_str());
        writer_vti->SetInputData(LIC);
        writer_vti->Write();

        std::string path_png = getDataPath("/data/test" + std::to_string(i) + ".png");
        vtkNew<vtkPNGWriter> writer_png;
        writer_png->SetFileName(path_png.c_str());
        writer_png->SetInputData(LIC);
        writer_png->Write();
    }
    /*
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
    */
    return EXIT_SUCCESS;
}