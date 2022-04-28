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

## Pressure Isolines: file pressure_iso.cpp
(Frawa and Annika)
Code displays a few isolines and has slider that can interactively display specific isolines. Data is weirdly rescaled rn but we have to ask the Prof. about it before we can do any better.
Needed files: ```./data/pressure_slice.vti```