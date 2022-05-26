//
// Created by peter on 07.05.2022.
//

#ifndef VIS2021_PLANEWIDGETINTERACTION_H
#define VIS2021_PLANEWIDGETINTERACTION_H

#include "utils.h"
#include "vtkCommand.h"
#include "vtkImageData.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPlaneWidget.h"
#include "vtkObject.h"

class PlaneWidgetInteraction : public vtkCommand {
public:
vtkTypeMacro(PlaneWidgetInteraction, vtkCommand);

    static PlaneWidgetInteraction *New() {
        return new PlaneWidgetInteraction;
    }


    void Execute(vtkObject *caller, unsigned long eventId, void *callData) {
        this->SetWidgetNormal();
        //Make sure plane is always in outline space of data (using pressure)
        UpdateWidget();

        updated = true;
        if (eventId == vtkCommand::EndInteractionEvent) {
            app->OnLeftButtonUp();
        }
    }

    void UpdateWidget() {
        double center[3];
        dataSpace->GetCenter(center);
        double bounds[6];
        dataSpace->GetBounds(bounds);

        double planeCenter[3];
        widget->GetCenter(planeCenter);
        planeCenter[0] = fmax(fmin(planeCenter[0], bounds[1]), bounds[0]);
        planeCenter[1] = fmax(fmin(planeCenter[1], bounds[3]), bounds[2]);
        planeCenter[2] = fmax(fmin(planeCenter[2], bounds[5]), bounds[4]);

        widget->SetCenter(planeCenter);

        if (zAxis) {
            widget->SetOrigin(bounds[0], bounds[2], planeCenter[2]);
            widget->SetPoint1(bounds[0], bounds[3], planeCenter[2]);
            widget->SetPoint2(bounds[1], bounds[2], planeCenter[2]);
        }
        if (xAxis) {
            widget->SetOrigin(planeCenter[0], bounds[2], bounds[4]);
            widget->SetPoint1(planeCenter[0], bounds[2], bounds[5]);
            widget->SetPoint2(planeCenter[0], bounds[3], bounds[4]);
        }
        if (yAxis) {
            widget->SetOrigin(bounds[0], planeCenter[1], bounds[4]);
            widget->SetPoint1(bounds[0], planeCenter[1], bounds[5]);
            widget->SetPoint2(bounds[1], planeCenter[1], bounds[4]);
        }
    }

    void SetWidgetNormal() {

        if (xAxis) {
            widget->SetNormal(1, 0, 0);
        }
        if (yAxis) {
            widget->SetNormal(0, 1, 0);
        }
        if (zAxis) {
            widget->SetNormal(0, 0, 1);
        }

    }

    void ToggleOrientation() {

        if (zAxis) {
            zAxis = false;
            xAxis = true;
        } else if (xAxis) {
            xAxis = false;
            yAxis = true;
        } else {
            zAxis = true;
            yAxis = false;
        }


        UpdateWidget();
    }

    bool zAxis = true;
    bool xAxis = false;
    bool yAxis = false;

    bool updated = false;
    vtkSmartPointer<vtkImageData> dataSpace;
    vtkSmartPointer<vtkPlaneWidget> widget;
    vtkSmartPointer<vtkInteractorStyle> app;
};

#endif //VIS2021_PLANEWIDGETINTERACTION_H
