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
        vtkPlaneWidget *widget = static_cast<vtkPlaneWidget *>(caller);
        widget->SetNormal(0, 0, 1);

        //Make sure plane is always in outline space of data (using pressure)
        double center[3];
        dataSpace->GetCenter(center);
        double bounds[6];
        dataSpace->GetBounds(bounds);

        double planeCenter[3];
        widget->GetCenter(planeCenter);
        planeCenter[0] = center[0];
        planeCenter[1] = center[1];

        planeCenter[2] = fmax(fmin(planeCenter[2], bounds[5]), bounds[4]);

        widget->SetCenter(planeCenter);

        widget->SetOrigin(bounds[0], bounds[2], planeCenter[2]);
        widget->SetPoint1(bounds[0], bounds[3], planeCenter[2]);
        widget->SetPoint2(bounds[1], bounds[2], planeCenter[2]);

        updated = true;
        if (eventId == vtkCommand::EndInteractionEvent) {
            app->OnLeftButtonUp();
        }
    }

    bool updated = false;
    vtkSmartPointer<vtkImageData> dataSpace;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> app;
};
#endif //VIS2021_PLANEWIDGETINTERACTION_H
