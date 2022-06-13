# meteo_cab


## Set up:
(Peter)
Place this repository as a subproject of the scivis repository.

Add ```add_subdirectory(meteo_cab)``` to the CmakeLists of the parent directory of this 
(just to have full build support)

Place our dataset at ```./data```

Check if build works -> Should show the same thing as prog03_volren

## Terrain: file terrainVis.cpp
(Annika)
The terrain is rendered by reading in a .obj mesh and the height is additionally read in via a .txt file to assign the colors (maybe this could be changed later though, this is just what currently works). The method is quite slow for some reason as of now.
The files can be found on the polybox folder ```https://polybox.ethz.ch/index.php/s/rW5ULLZJKK7SMdK``` and can be placed into the data folder.


## Cloud rendering:
(Peter)

cloudyness.mhd -> Sum of CLI and CLW properties. Rendered through raycasting.

File can be replaced any time with a higher resolution sampling.
Load Cloudyness.pvsm -> Set Resample dimensions -> Export -> enjoy

## Pressure Slices: file slice_pressure.cpp
(Frawa and Annika)

Code displays a few isolines, arrow up and down chooses different slice and left and right different isolines values. Pressing "h" switches between horizontal and vertical slices. Horizontal coloring is done locally to see more change and has adaptive legend. Vertical coloring on the other hand is done globally because it is not very informative anyways (and thus also has no isolines, but if you'd like we could change that)

Needed files: ```./data/pressure_slice.vti```



## Iterate through Slices of ImageData:

```
    vtkNew<vtkImageReslice> slicer;
    double spacing[3];
    data->GetSpacing(spacing);
    double origin[3];
    double cur_origin[3];
    data->GetOrigin(origin);
    data->GetOrigin(cur_origin);
    slicer->SetInputData(data);
    int extent[3];
    data->GetExtent(extent);
    slicer->SetOutputExtent(0, extent[0] - 1, 0, extent[1] - 1, 0, 0);
    for (int i = 0; i < extent[2]; ++i) {
        cur_origin[2] = origin[2] + i * spacing[2];
        slicer->SetOutputOrigin(cur_origin);
        slicer->Update();
        vtkSmartPointer<vtkImageData> slice = slicer->GetOutput();

        //do your stuff
    }
```
