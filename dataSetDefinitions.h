//
// Created by peter on 13.04.2022.
//

#ifndef VIS2021_DATASETDEFINITIONS_H
#define VIS2021_DATASETDEFINITIONS_H


#include <string>
#include <algorithm>

#include <cstdlib>

static std::string getDataPath(const std::string &localDataPath) {
    auto basepath = std::string(__FILE__);
    std::replace(basepath.begin(), basepath.end(), '\\', '/');
    auto lastSlash = basepath.find_last_of('/');
    basepath = basepath.substr(0, lastSlash);
    auto fullpath = basepath + localDataPath;
    return fullpath;
}


#endif //VIS2021_DATASETDEFINITIONS_H
