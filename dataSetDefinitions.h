//
// Created by peter on 13.04.2022.
//

#ifndef VIS2021_DATASETDEFINITIONS_H
#define VIS2021_DATASETDEFINITIONS_H


#include <string>

using namespace std;
enum DataSets {
    CLOUD_ICE,
    CLOUD_WATER,
    PRESSURE,
    RAIN_MIXTURE,
    WIND_VERTICAL,
    WIND_LONGITUDIAL,
    WIND_LATERAL
};

static string getDatasetPath(DataSets dataSets) {
    string base = "../../../data/";
    switch (dataSets) {
        case DataSets::CLOUD_ICE:
            base += "cli";
            break;
        case DataSets::CLOUD_WATER:
            base += "clw";
            break;
        case DataSets::PRESSURE:
            base += "pres";
            break;
        case DataSets::RAIN_MIXTURE:
            base += "qr";
            break;
        case DataSets::WIND_VERTICAL:
            base += "ua";
            break;
        case DataSets::WIND_LATERAL:
            base += "va";
            break;
        case DataSets::WIND_LONGITUDIAL:
            base += "wa";
            break;
    }
    base += "/";
    return base;
}


#endif //VIS2021_DATASETDEFINITIONS_H
