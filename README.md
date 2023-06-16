Instructions

#Make sure to have .flat.root files in your mySamples directory
#To make them, go over to https://github.com/maxdubnowski/BuildEventGenerators.git and follow this Repo
#Transfer the .flat.root files to the mySamples directory

#Run the Flat tree analyzer which loops through the .flat files,  makes all of the plots, and converts to .root file for each of the event generators, 

root -l script_LoopGenerators.cxx


#Plot the 1D histograms

root -l GeneratorOverlay.cpp


#Plot the 2D histograms

root -l GeneratorOverlay2D.cpp
