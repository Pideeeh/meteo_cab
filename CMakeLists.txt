cmake_minimum_required(VERSION 3.8...3.12 FATAL_ERROR)

PROJECT(meteocab)

find_package(VTK REQUIRED)

add_executable(${PROJECT_NAME} MACOSX_BUNDLE
        ${PROJECT_NAME}.cpp
        vtkEasyTransfer.hpp
        dataSetDefinitions.h)

target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})

add_executable(Terrain_Vis terrainVis.cpp)
target_link_libraries(Terrain_Vis PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS Terrain_Vis MODULES ${VTK_LIBRARIES})


add_executable(Pressure_Iso pressure_iso.cpp)
target_link_libraries(Pressure_Iso PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS Pressure_Iso MODULES ${VTK_LIBRARIES})

add_executable(Slicer_Pressure slicer_pressure.cpp)
target_link_libraries(Slicer_Pressure PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS Slicer_Pressure MODULES ${VTK_LIBRARIES})


add_executable(Wind_2d wind_2d.cpp)
target_link_libraries(Wind_2d PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS Wind_2d MODULES ${VTK_LIBRARIES})

add_executable(Resampling resampling.cpp)
target_link_libraries(Resampling PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS Resampling MODULES ${VTK_LIBRARIES})


add_executable(App MeteoAppLauncher.cpp)
target_link_libraries(App PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS App MODULES ${VTK_LIBRARIES})


add_executable(LICbuilder LIC_builder.cpp)
target_link_libraries(LICbuilder PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS LICbuilder MODULES ${VTK_LIBRARIES})