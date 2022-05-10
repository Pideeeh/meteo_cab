#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkSmartPointer.h"
//new
#include "vtkXMLImageDataReader.h" from 'C:\Users\annik\Documents\fs22\scientific visualization\vis22_prog\VTK-9.1.0\build\IO\XML\IOXML.dir\Release\vtkXMLImageDataReader.obj'
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkWarpScalar.h>
#include <iostream>
#include <fstream>
#include <vtkXMLImageDataWriter.h>

using namespace std;

//gives slice indices and coeff for upper
void get_surrounding_slices(double height, int& upper, int& lower, double& coeff) {
	ifstream myfile;
	myfile.open("C:/Users/annik/Documents/fs22/scientific visualization/vis22_prog/build/meteo_cab/Resampling.dir/Debug/sampling_height.txt");
	string str;
	if (height <= 106.81226) {//lowest part
		upper = 0;
		lower = 0;
		coeff = 1;
		return;
	}
	if (height >= 20808.89648) {
		upper = 149;
		lower = 149;
		coeff = 1;
		return;
	}
	upper = 0;
	lower = 0;
	double upper_val, lower_val;
	for (int i = 0; i < 150; i++) {
		getline(myfile, str);
		double sample_height = stod(str) * 0.0001;
		if (sample_height < height) {
			lower = i;
			lower_val = sample_height;
			break;
		}
		upper = i;
		upper_val = sample_height;
	}
	coeff = (height - lower_val) / (upper_val - lower_val);
}

int main()
{
	vtkNew<vtkXMLImageDataReader> reader;
	reader->SetFileName("C:/Users/annik/Documents/fs22/scientific visualization/vis22_prog/meteo_cab/data/pres/pres_10.vti_scaled.vti");
	reader->Update();
	vtkSmartPointer<vtkImageData> dat = vtkSmartPointer<vtkImageData>::New();
	dat = reader->GetOutput();
	int* dims = dat->GetDimensions();
	vtkSmartPointer<vtkPointData> pts = dat->GetPointData();
	vtkSmartPointer<vtkAbstractArray> arr = pts->GetAbstractArray("pres");
	vtkSmartPointer<vtkFloatArray> float_array = vtkFloatArray::SafeDownCast(arr);
	cout << "z dim is " << dims[2] << endl;
	cout << "y dim is " << dims[1] << endl;
	cout << "x dim is " << dims[0] << endl;

	int resolution = 100;
	int new_dims[3];
	new_dims[0] = dims[0];
	new_dims[1] = dims[1];
	new_dims[2] = resolution;
	
	vtkSmartPointer<vtkImageData> regular_data = vtkSmartPointer<vtkImageData>::New();
	regular_data->SetDimensions(new_dims);
	auto bd=regular_data->GetBounds();
	bd[0] = dat->GetBounds()[0];
	bd[1] = dat->GetBounds()[1];
	bd[2] = dat->GetBounds()[2];
	bd[3] = dat->GetBounds()[3];
	bd[4] = dat->GetBounds()[4];
	bd[5] = dat->GetBounds()[5];
	regular_data->AllocateScalars(VTK_DOUBLE, 1);

	//vtkSmartPointer<vtkPointData> reg_pts = regular_data->GetPointData();
	//vtkSmartPointer<vtkAbstractArray> reg_arr = reg_pts->GetAbstractArray("pres");
	//vtkSmartPointer<vtkFloatArray> reg_float_array = vtkFloatArray::SafeDownCast(reg_arr);

	double upper_height = 20808*0.0001;
	double lower_height = 107* 0.0001;

	for (int z = 0; z < new_dims[2]; z++) {
		cout << z << endl;
		ifstream myfile;
		myfile.open("C:/Users/annik/Documents/fs22/scientific visualization/vis22_prog/build/meteo_cab/Resampling.dir/Debug/height_values.txt");
		string str;
		for (int y = 0; y < new_dims[1]; y++) {
			cout << y << endl;
			for (int x = 0; x < new_dims[0]; x++) {
				getline(myfile, str);
				double terrain_height = stod(str) * 0.0001;
				double regular_sample_in_data_space = z*(upper_height-lower_height)/resolution + lower_height + terrain_height;
				int upper, lower;
				double coeff;
				get_surrounding_slices(regular_sample_in_data_space, upper, lower, coeff);
				int idx_dat = upper * dims[0] * dims[1] + y * dims[0] + x;
				double upper_val= static_cast<double*>(float_array->GetTuple(idx_dat))[0];
				idx_dat = lower * dims[0] * dims[1] + y * dims[0] + x;
				double lower_val = static_cast<double*>(float_array->GetTuple(idx_dat))[0];
				double new_dat_val = coeff * upper_val + (1.0 - coeff) * lower_val;
				double* v = static_cast<double*>(regular_data->GetScalarPointer(x, y, z));
				v[0]=new_dat_val;
			}
		}
	}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("data_regular.vti");
	writer->SetInputData(regular_data);
	writer->Write();

	//for (int z = 0; z < dims[2]; z++) {
	//	ifstream myfile;
	//	myfile.open("height_values.txt");
	//	string str;
	//	getline(myfile, str);
	//	double terrain_height=stod(str) * 0.0001;
	//	for (int y = 0; y < dims[1]; y++) {
	//		for (int x = 0; x < dims[0]; x++) {
	//			int idx = z*dims[0]*dims[1]+y * dims[0] + x;
	//			auto v = static_cast<double*>(float_array->GetTuple(idx));
	//			//myfile << v[0] << endl;
	//			double actual_given_sample_height=v[0]+terrain_height;
	//		}
	//	}
	//}
	return 0;
}