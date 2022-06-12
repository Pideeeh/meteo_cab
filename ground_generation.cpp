#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/writeOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/adjacency_list.h>
#include <igl/barycenter.h>
#include <igl/massmatrix.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/invert_diag.h>
#include <igl/arap.h>
#include <igl/per_vertex_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <cmath>
#include <igl/Timer.h>
//#include <array>

using namespace Eigen;
using namespace std;


//code to create mesh of ground 
void planegen(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int elementsx, int elementsy) {
    ifstream myfile("C:/Users/annik/Documents/fs22/scientific visualization/vis22_prog/meteo_cab/data/height_values.txt");
    if (!myfile.is_open())
    {
        cout << "error opening file" << endl;
        return;
    }
    else {
        cout << "file openend" << endl;
    }
    Eigen::VectorXd rangex = Eigen::VectorXd::LinSpaced(elementsx, 4.5, 14.496);
    Eigen::VectorXd rangey = Eigen::VectorXd::LinSpaced(elementsy, 47.5, 54.4975);
    int m = elementsx * elementsy; //number of vertices
    V = Eigen::MatrixXd::Zero(m, 3);
    int faces = 2 * (elementsx - 1) * (elementsy - 1);
    F = Eigen::MatrixXi::Zero(faces, 3);
    int cf = 0;
    int p;
    for (int y = 0; y < elementsy; y++) {
        for (int x = 0; x < elementsx; x++) {
            p = y * elementsx + x;
            V(p, 0) = rangex(x);
            V(p, 1) = rangey(y);
            string str;
            getline(myfile, str);
            V(p, 2) = stod(str) * 0.0001;
            if (x < elementsx - 1 && y < elementsy - 1) {//can create face from it
                F(cf, 2) = p;
                F(cf, 1) = p + elementsx;
                F(cf, 0) = p + 1;
                cf++;
            }
            if (x > 0 && y > 0) {
                F(cf, 2) = p;
                F(cf, 1) = p - elementsx;
                F(cf, 0) = p - 1;
                cf++;
            }
        }
    }
    myfile.close();
}

int main() {
    MatrixXd V;
    MatrixXi F;
    planegen(V, F, 1429, 1556);
    igl::writeOBJ("mesh.obj", V, F);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);//paint mesh. U is V, F faces
    viewer.launch();
    return 0;
}
