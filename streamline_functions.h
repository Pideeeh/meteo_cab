//
// Created by peter on 16.05.2022.
//

#ifndef VIS2021_STREAMLINE_FUNCTIONS_H
#define VIS2021_STREAMLINE_FUNCTIONS_H

#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkNamedColors.h"
#include "vtkLine.h"


//FROM ANNIKA AND FRAWA TO DAVIDE
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


//FROM ANNIKA AND FRAWA TO DAVIDE
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


void compute_one_streamline(vtkSmartPointer<vtkImageData> v, double *point, vtkSmartPointer<vtkPolyData> result) {

    //// Allocate one seed
    //vtkNew<vtkPolyData> seeds_LIC;
    //vtkNew<vtkPoints> one_seed;
    //one_seed->InsertNextPoint(point);
    //seeds_LIC->SetPoints(one_seed);

    //// Compute streamline
    //vtkNew<vtkStreamTracer> tracer_LIC;
    //tracer_LIC->SetInputData(v);
    //tracer_LIC->SetSourceData(seeds_LIC);
    //tracer_LIC->SetMaximumPropagation(1);
    //tracer_LIC->SetIntegratorTypeToRungeKutta45();
    //tracer_LIC->SetInitialIntegrationStep(0.01);
    //tracer_LIC->SetIntegrationDirectionToForward(); // TODO: should be to both for LIC texture
    //tracer_LIC->Update();

    //// return streamline points
    //vtkSmartPointer<vtkPoints> points = tracer_LIC->GetOutput()->GetPoints();
    //result->SetPoints(points);

    vtkNew<vtkPoints> line_points;
    int n_steps = 1000;
    double s = 0.01;

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


// Turn result of streamtraces (vtkPoints) into lines (which can be visualized with PolyDataMapper).
void from_points_to_line(
        vtkSmartPointer<vtkPoints> line_points,
        vtkSmartPointer<vtkPoints> all_points,
        vtkSmartPointer<vtkCellArray> all_lines,
        vtkSmartPointer< vtkUnsignedCharArray> all_colors
) {
    int n_points = line_points->GetNumberOfPoints();
    int offset = all_points->GetNumberOfPoints();

    vtkNew<vtkNamedColors> colors;

    for (int i = 0; i < n_points-1; i++) {
        // create a segment between each set of point
        double *point = line_points->GetPoint(i);
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, offset + i);
        line->GetPointIds()->SetId(1, offset + i + 1);

        // update set of all points, colors, and lines
        all_points->InsertNextPoint(point);
        all_lines->InsertNextCell(line);
        all_colors->InsertNextTuple(colors->GetColor3d("Blue").GetData());
    }

    // Add the last point to the array
    double* point = line_points->GetPoint(n_points-1);
    all_points->InsertNextPoint(point);
}

#endif //VIS2021_STREAMLINE_FUNCTIONS_H
