//
// Created by peter on 05.05.2022.
//

#include "MeteoCabApp.h"

int main(int, char *[]) {

    MeteoCabApp app{};

    cout << "Press '1' to toggle pressure and move the slice vertically with your mouse" << endl;
    cout << "Press '2' to toggle the terrain" << endl;
    cout << "Press '3' to toggle the cloud rendering" << endl;
    cout << "Press '4' to toggle the wind and move the slice vertically with your mouse" << endl;
    cout << "Press '5' to toggle the divergence and move the slice vertically with your mouse" << endl;

    app.Launch();

}