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
#include "vtkXMLImageDataReader.h"
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkWarpScalar.h>
#include <iostream>
#include <fstream>
#include <vtkXMLImageDataWriter.h>

#include "dataSetDefinitions.h"

using namespace std;

//gives slice indices and coeff for upper
void get_surrounding_slices(double height, int& upper, int& lower, double& coeff, std::vector<double>& sample_height) {
	if (height <= sample_height[sample_height.size()-1]) {//lowest part
		upper = sample_height.size() - 1;
		lower = sample_height.size() - 1;
		coeff = 1;
		return;
	}
	if (height >= sample_height[0]) {
		upper = 0;
		lower = 0;
		coeff = 1;
		return;
	}
	upper = 0;
	lower = 0;
	double upper_val, lower_val;
	for (int i = 0; i < sample_height.size(); i++) {
		if (sample_height[i] < height) {
			lower = i;
			lower_val = sample_height[i];
			break;
		}
		upper = i;
		upper_val = sample_height[i];
	}
	coeff = (height - lower_val) / (upper_val - lower_val);
}

int main()
{
	vtkSmartPointer<vtkXMLImageDataReader> reader=vtkSmartPointer<vtkXMLImageDataReader>::New();
	//----------------------------
	// this is dataset specific
	//---------------------------------
	reader->SetFileName(getDataPath("/data/pres/pres_10.vti_scaled.vti").c_str());//adapt
	//--------------------------
	reader->Update();
	vtkSmartPointer<vtkImageData> dat = vtkSmartPointer<vtkImageData>::New();
	dat = reader->GetOutput();
	int* dims = dat->GetDimensions();
	vtkSmartPointer<vtkPointData> pts = dat->GetPointData();
	//----------------------------
	// this is dataset specific
	//---------------------------------
	vtkSmartPointer<vtkAbstractArray> arr = pts->GetAbstractArray("pres");//adapt
	//----------------------------------
	vtkSmartPointer<vtkFloatArray> float_array = vtkFloatArray::SafeDownCast(arr);

	int resolution = 200; //adapt: choose how many samples along the z axis we want to have.
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
	ifstream myfile;
	myfile.open(getDataPath("/data/height.txt").c_str());//For wa, adapt to height2.txt. For all others, adapt to height.txt
	string str;
	std::vector<double> height;
	int num_height_samples = 150; //For wa, adapt to 151. For all others, adapt to 150. 
	for (int i = 0; i < num_height_samples; i++) {
		getline(myfile, str);
		height.push_back(stod(str) * 0.0001);
	}
	double upper_height = 20808*0.0001;
	double lower_height = 107*0.00001;

	regular_data->SetSpacing(dat->GetSpacing()[0], dat->GetSpacing()[1], -(upper_height - lower_height) / resolution);
	regular_data->SetOrigin(dat->GetOrigin()[0], dat->GetOrigin()[1], dat->GetOrigin()[2]);

	
	myfile.open(getDataPath("/data/height_values.txt").c_str());
	for (int y = 0; y < new_dims[1]; y++) {
		for (int x = 0; x < new_dims[0]; x++) {
			getline(myfile, str);
			double terrain_height = stod(str) * 0.0001;
			for (int z = 0;z<new_dims[2]; z++) {
				double regular_sample_in_data_space = (z) * (lower_height - upper_height) / resolution + upper_height - terrain_height;
				int upper, lower;
				double coeff;
				get_surrounding_slices(regular_sample_in_data_space, upper, lower, coeff, height);
				int idx_dat = upper * dims[0] * dims[1] + y * dims[0] + x;
				double upper_val = static_cast<double*>(float_array->GetTuple(idx_dat))[0];
				idx_dat = lower * dims[0] * dims[1] + y * dims[0] + x;
				double lower_val = static_cast<double*>(float_array->GetTuple(idx_dat))[0];
				double new_dat_val = coeff * upper_val + (1.0 - coeff) * lower_val;
				double* v = static_cast<double*>(regular_data->GetScalarPointer(x,y, z));
				v[0] = new_dat_val;
			}
		}
	}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(getDataPath("/data/regular_data.vti").c_str());
	writer->SetInputData(regular_data);
	writer->Write();
	return 0;
}