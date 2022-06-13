# meteo_cab


Our cpp-Codes as well as some paraview state files can be found in this github repository. Our main App can be Launched via MeteoAppLauncher.cpp and it contains most of our components and calculations. The resampling of the data however we did beforehand via resampling.cpp, the high-res terrain mesh was generated externally using libigl in the ground_generation.cpp file (after reading out the z_ifc data one by one and converting it to a text document) and then downsampled using Meshlab. The relevant results can be found under ```https://polybox.ethz.ch/index.php/s/rW5ULLZJKK7SMdK```. The code wind_2d.cpp is relevant for doing the wind computations and also includes our LIC texture generation, which we used to save the LIC textures beforehand and just loaded them in the App to allow for better runtime.

###Setup
Place this repository as a subproject of the scivis repository.

Add ```add_subdirectory(meteo_cab)``` to the CmakeLists of the parent directory of this 
(just to have full build support)

Place our dataset at ```./data```
