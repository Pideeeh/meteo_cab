# meteo_cab


## Set up:
Place this repository as a subproject of the scivis repository.

Add ```add_subdirectory(meteo_cab)``` to the CmakeLists of the parent directory of this 
(just to have full build support)

Place our dataset at ```./data```

Check if build works -> Should show the same thing as prog03_volren

## Terrain: file terrainVis.cpp
The terrain is rendered by reading in a .obj mesh and the height is additionally read in via a .txt file to assign the colors (maybe this could be changed later though, this is just what currently works). The method is quite slow for some reason as of now.
The files can be found on the polybox folder ```https://polybox.ethz.ch/index.php/s/rW5ULLZJKK7SMdK``` and can be placed into the data folder.
ATTENTION: the .cpp file still contains *my* path to the data as the one from Peter did not work for me, with the environment variable this can be adapted so we can all use the same code:) But for now, you just have to change the two paths.


## Cloud rendering:

cloudyness.mhd -> Sum of CLI and CLW properties. Rendered through raycasting.

File can be replaced any time with a higher resolution sampling.
Load Cloudyness.pvsm -> Set Resample dimensions -> Export -> enjoy