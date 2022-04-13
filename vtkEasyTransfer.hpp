#ifndef INCLUDE_ONCE_VTK_EASY_TRANSFER
#define INCLUDE_ONCE_VTK_EASY_TRANSFER

#include <vtkObject.h>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkImageMapper.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkActor2D.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <algorithm>

// Very simplistic transfer function editor for VTK.
class vtkEasyTransfer : public vtkObject
{
public:
	friend class vtkEasyTransferInteractorStyle;

	static vtkEasyTransfer* New();

	// Constructor.
	vtkEasyTransfer();

	// Sets a default color map.
	void SetColormapHeat()
	{
		mColorTransferFunction->RemoveAllPoints();
		mOpacityTransferFunction->RemoveAllPoints();
		int numPoints = 256;
		double minValue = 0, maxValue = 255;
		for (int i = 0; i < numPoints; ++i)
		{
			double t = i / (double)(numPoints - 1);
			double x = minValue + t * (maxValue - minValue);
			double r = std::min(std::max(0.0, ((-1.0411*t + 0.4149)*t + 0.0162)*t + 0.9949), 1.0);
			double g = std::min(std::max(0.0, ((1.7442*t + -2.7388)*t + 0.1454)*t + 0.9971), 1.0);
			double b = std::min(std::max(0.0, ((0.9397*t + -0.2565)*t + -1.5570)*t + 0.9161), 1.0);
			mColorTransferFunction->AddRGBPoint(x, r, g, b);
			mOpacityTransferFunction->AddPoint(x, t);
		}
		RefreshImage();
	}

	// Sets a default color map.
	void SetColormapBlue2Red()
	{
		mColorTransferFunction->RemoveAllPoints();
		mOpacityTransferFunction->RemoveAllPoints();
		int numPoints = 256;
		double minValue = 0, maxValue = 255;
		for (int i = 0; i < numPoints; ++i)
		{
			double t = i / (double)(numPoints - 1);
			double x = minValue + t * (maxValue - minValue);
			double r = std::min(std::max(0.0, (((5.0048*t + -8.0915)*t + 1.1657)*t + 1.4380)*t + 0.6639), 1.0);
			double g = std::min(std::max(0.0, (((7.4158*t + -15.9415)*t + 7.4696)*t + 1.2767)*t + -0.0013), 1.0);
			double b = std::min(std::max(0.0, (((6.1246*t + -16.2287)*t + 11.9910)*t + -1.4886)*t + 0.1685), 1.0);
			mColorTransferFunction->AddRGBPoint(x, r, g, b);
			mOpacityTransferFunction->AddPoint(x, t);
		}
		RefreshImage();
	}

	// Sets the color map range
	void SetColorRange(double minValue, double maxValue)
	{
		double range[2] = { mColorTransferFunction->GetRange()[0], mColorTransferFunction->GetRange()[1] };
		double* dataPtr = mColorTransferFunction->GetDataPointer();
		int size = mColorTransferFunction->GetSize();
		std::vector<double> data(dataPtr, dataPtr + (4 * size));
		mColorTransferFunction->RemoveAllPoints();
		for (int i = 0; i < size; ++i)
		{
			double told = (data[i * 4] - range[0]) / (range[1] - range[0]);
			double tnew = told * (maxValue - minValue) + minValue;
			mColorTransferFunction->AddRGBPoint(tnew, data[i * 4 + 1], data[i * 4 + 2], data[i * 4 + 3]);
		}
	}

	// Sets the opacity map range
	void SetOpacityRange(double minValue, double maxValue)
	{
		double range[2] = { mOpacityTransferFunction->GetRange()[0], mOpacityTransferFunction->GetRange()[1] };
		double* dataPtr = mOpacityTransferFunction->GetDataPointer();
		int size = mOpacityTransferFunction->GetSize();
		std::vector<double> data(dataPtr, dataPtr + (2 * size));
		mOpacityTransferFunction->RemoveAllPoints();
		for (int i = 0; i < size; ++i)
		{
			double told = (data[i * 2] - range[0]) / (range[1] - range[0]);
			double tnew = told * (maxValue - minValue) + minValue;
			mOpacityTransferFunction->AddPoint(tnew, data[i * 2 + 1]);
		}
	}

	// Gets the color transfer function.
	vtkSmartPointer<vtkColorTransferFunction> GetColorTransferFunction() { return mColorTransferFunction; }
	// Gets the opacity transfer function.
	vtkSmartPointer<vtkPiecewiseFunction> GetOpacityTransferFunction() { return mOpacityTransferFunction; }
	// Gets the renderer.
	vtkSmartPointer<vtkRenderer> GetRenderer() { return mRenderer; }
	// Gets the interactor for user interactions.
	vtkSmartPointer<vtkInteractorStyle> GetInteractorStyle() { return mInteractorStyle; }

	// Recomputes the color map images. Call this function, whenever the viewport has changed! (only updates the viewport bounds, if this is already assigned to a vtkRenderWindow)
	void RefreshImage()
	{
		if (mRenderer->GetRenderWindow())
		{
			int* size = mRenderer->GetRenderWindow()->GetSize();
			double* viewport = mRenderer->GetViewport();
			mCanvasWidth = size[0] * (viewport[2] - viewport[0]);
			mCanvasHeight = size[1] * (viewport[3] - viewport[1]) / 4;
			mDrawing->SetExtent(0, size[0] * (viewport[2] - viewport[0]) - 1, 0, size[1] * (viewport[3] - viewport[1]) - 1, 0, 1);
		}

		mDrawing->SetDrawColor(1, 1, 1);
		mDrawing->FillBox(0, mCanvasWidth, 0, 4 * mCanvasHeight);

		// draw color maps
		{
			double range[2] = { mColorTransferFunction->GetRange()[0], mColorTransferFunction->GetRange()[1] };
			double* dataPtr = mColorTransferFunction->GetDataPointer();
			int size = mColorTransferFunction->GetSize();
			std::vector<double> data(dataPtr, dataPtr + (4 * size));

			for (int i = 0; i < size; ++i)
			{
				int lower = (data[4 * i] - range[0]) / (range[1] - range[0]) * (mCanvasWidth - 1);
				int upper = i == size - 1 ? mCanvasWidth : ((data[4 * i + 4] - range[0]) / (range[1] - range[0]) * (mCanvasWidth - 1));
				int valueR = data[4 * i + 1] * mCanvasHeight;
				int valueG = data[4 * i + 2] * mCanvasHeight;
				int valueB = data[4 * i + 3] * mCanvasHeight;

				mDrawing->SetDrawColor(1.0, 0.1, 0.1);
				mDrawing->FillBox(lower, upper, 3 * mCanvasHeight, 3 * mCanvasHeight + valueR);
				mDrawing->SetDrawColor(0.1, 1.0, 0.1);
				mDrawing->FillBox(lower, upper, 2 * mCanvasHeight, 2 * mCanvasHeight + valueG);
				mDrawing->SetDrawColor(0.1, 0.1, 1.0);
				mDrawing->FillBox(lower, upper, 1 * mCanvasHeight, 1 * mCanvasHeight + valueB);
			}
		}
		// draw opacity map
		{
			double range[2] = { mOpacityTransferFunction->GetRange()[0], mOpacityTransferFunction->GetRange()[1] };
			double* dataPtr = mOpacityTransferFunction->GetDataPointer();
			int size = mOpacityTransferFunction->GetSize();
			std::vector<double> data(dataPtr, dataPtr + (2 * size));

			for (int i = 0; i < size; ++i)
			{
				int lower = (data[2 * i] - range[0]) / (range[1] - range[0]) * (mCanvasWidth - 1);
				int upper = i == size - 1 ? mCanvasWidth : ((data[2 * i + 2] - range[0]) / (range[1] - range[0]) * (mCanvasWidth - 1));
				int valueO = data[2 * i + 1] * mCanvasHeight;

				mDrawing->SetDrawColor(0.3, 0.3, 0.3);
				mDrawing->FillBox(lower, upper, 0 * mCanvasHeight, 0 * mCanvasHeight + valueO);
			}
		}
	}

private:
	vtkEasyTransfer(const vtkEasyTransfer&); // Not implemented.
	void operator=(const vtkEasyTransfer&); // Not implemented.

	vtkSmartPointer<vtkColorTransferFunction> mColorTransferFunction;
	vtkSmartPointer<vtkPiecewiseFunction> mOpacityTransferFunction;
	vtkSmartPointer<vtkRenderer> mRenderer;
	vtkSmartPointer<vtkImageCanvasSource2D> mDrawing;
	vtkSmartPointer<vtkEasyTransferInteractorStyle> mInteractorStyle;

	int mCanvasWidth;
	int mCanvasHeight;

};
vtkStandardNewMacro(vtkEasyTransfer);

// Interactor style
class vtkEasyTransferInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static vtkEasyTransferInteractorStyle* New();
	vtkEasyTransferInteractorStyle() : mEasyTransfer(NULL), mMouseDown(false), mLastX(-1), mLastPlot(-1) {}
	void SetEasyTransfer(vtkEasyTransfer* easyTransfer) { mEasyTransfer = easyTransfer; }
	virtual void OnLeftButtonDown() override {
		int* clickPos = this->GetInteractor()->GetEventPosition();
		mMouseDown = true;
		int* size = mEasyTransfer->mRenderer->GetRenderWindow()->GetSize();
		double windowRelative[2] = { clickPos[0] / (double)size[0], clickPos[1] / (double)size[1] };
		double* viewport = mEasyTransfer->mRenderer->GetViewport();
		double viewportRelative[2] = {
			(windowRelative[0] - viewport[0]) / (viewport[2] - viewport[0]),
			(windowRelative[1] - viewport[1]) / (viewport[3] - viewport[1])
		};
		if (0 <= viewportRelative[0] && viewportRelative[0] <= 1 && 0 <= viewportRelative[1] && viewportRelative[1] <= 1)
		{
			int x = viewportRelative[0] * mEasyTransfer->mCanvasWidth;
			mLastX = -1;
			int yfull = viewportRelative[1] * (mEasyTransfer->mCanvasHeight * 4);
			int plot = yfull / mEasyTransfer->mCanvasHeight;
			mLastPlot = plot;
			processPixel(clickPos);
		}
		else mLastPlot = -1;
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}
	virtual void OnLeftButtonUp() override {
		mMouseDown = false;
		vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	}
	virtual void OnMouseMove() override {
		if (mMouseDown)
		{
			int* clickPos = this->GetInteractor()->GetEventPosition();
			processPixel(clickPos);
		}
		vtkInteractorStyleTrackballCamera::OnMouseMove();
	}
private:
	bool mMouseDown;
	int mLastX;
	int mLastPlot;
	void processPixel(int* clickPos)
	{
		int* size = mEasyTransfer->mRenderer->GetRenderWindow()->GetSize();
		double windowRelative[2] = { clickPos[0] / (double)size[0], clickPos[1] / (double)size[1] };
		double* viewport = mEasyTransfer->mRenderer->GetViewport();
		double viewportRelative[2] = {
			(windowRelative[0] - viewport[0]) / (viewport[2] - viewport[0]),
			(windowRelative[1] - viewport[1]) / (viewport[3] - viewport[1])
		};
		if (0 <= viewportRelative[0] && viewportRelative[0] <= 1 && 0 <= viewportRelative[1] && viewportRelative[1] <= 1)
		{
			int x = viewportRelative[0] * mEasyTransfer->mCanvasWidth;
			int yfull = viewportRelative[1] * (mEasyTransfer->mCanvasHeight * 4);
			int plot = yfull / mEasyTransfer->mCanvasHeight;
			if (mLastPlot != plot) return;
			int y = yfull % mEasyTransfer->mCanvasHeight;
			// color plots: 3=R, 2=G, 1=B
			if (plot >= 1) {
				double range[2] = { mEasyTransfer->mColorTransferFunction->GetRange()[0], mEasyTransfer->mColorTransferFunction->GetRange()[1] };
				double* dataPtr = mEasyTransfer->mColorTransferFunction->GetDataPointer();
				int size = mEasyTransfer->mColorTransferFunction->GetSize();
				std::vector<double> data(dataPtr, dataPtr + (4 * size));
				double T = viewportRelative[0] * (range[1] - range[0]) + range[0];

				int L = 0, R = size;
				while (L < R) {
					int m = (int)((L + R) / 2);
					if (data[4 * m] < T)
						L = m + 1;
					else R = m;
				}
				L = std::max(L - 1, 0);
				data[4 * L + 4 - plot] = y / (double)mEasyTransfer->mCanvasHeight;
				if (mLastX == -1) mLastX = L;
				for (int index = mLastX; index <= L + 1; ++index)
					mEasyTransfer->mColorTransferFunction->AddRGBPoint(data[4 * index], data[4 * (plot == 3 ? L : index) + 1], data[4 * (plot == 2 ? L : index) + 2], data[4 * (plot == 1 ? L : index) + 3]);
				for (int index = L; index <= mLastX + 1; ++index)
					mEasyTransfer->mColorTransferFunction->AddRGBPoint(data[4 * index], data[4 * (plot == 3 ? L : index) + 1], data[4 * (plot == 2 ? L : index) + 2], data[4 * (plot == 1 ? L : index) + 3]);
				mLastX = L;
				mEasyTransfer->RefreshImage();
			}
			// opacity plot: 0
			else {
				double range[2] = { mEasyTransfer->mOpacityTransferFunction->GetRange()[0], mEasyTransfer->mOpacityTransferFunction->GetRange()[1] };
				double* dataPtr = mEasyTransfer->mOpacityTransferFunction->GetDataPointer();
				int size = mEasyTransfer->mOpacityTransferFunction->GetSize();
				std::vector<double> data(dataPtr, dataPtr + (2 * size));
				double T = viewportRelative[0] * (range[1] - range[0]) + range[0];

				int L = 0, R = size;
				while (L < R) {
					int m = (int)((L + R) / 2);
					if (data[2 * m] < T)
						L = m + 1;
					else R = m;
				}
				L = std::max(L - 1, 0);
				data[2 * L + 1] = y / (double)mEasyTransfer->mCanvasHeight;
				if (mLastX == -1) mLastX = L;
				for (int index = mLastX; index <= L + 1; ++index)
					mEasyTransfer->mOpacityTransferFunction->AddPoint(data[2 * index], data[2 * L + 1]);
				for (int index = L; index <= mLastX + 1; ++index)
					mEasyTransfer->mOpacityTransferFunction->AddPoint(data[2 * index], data[2 * L + 1]);
				mLastX = L;
				mEasyTransfer->RefreshImage();
			}
		}
	}
	vtkEasyTransfer* mEasyTransfer;
};
vtkStandardNewMacro(vtkEasyTransferInteractorStyle);

vtkEasyTransfer::vtkEasyTransfer() : mCanvasWidth(256), mCanvasHeight(100)
{
	mColorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	mOpacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();

	mRenderer = vtkSmartPointer<vtkRenderer>::New();
	mRenderer->ResetCamera();

	mDrawing = vtkSmartPointer<vtkImageCanvasSource2D>::New();
	mDrawing->SetNumberOfScalarComponents(3);
	mDrawing->SetScalarTypeToUnsignedChar();
	mDrawing->SetExtent(0, mCanvasWidth - 1, 0, 4 * mCanvasHeight - 1, 0, 1);

	vtkSmartPointer<vtkImageMapper> imageMapper = vtkSmartPointer<vtkImageMapper>::New();
	imageMapper->SetInputConnection(mDrawing->GetOutputPort());
	imageMapper->SetColorWindow(1);
	imageMapper->SetColorLevel(0);

	vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
	imageActor->SetMapper(imageMapper);
	mRenderer->AddActor2D(imageActor);

	mInteractorStyle = vtkSmartPointer<vtkEasyTransferInteractorStyle>::New();
	mInteractorStyle->SetEasyTransfer(this);
}

#endif