#!/bin/bash

# Creation of the prototype for the hippocampus:
../../utilities/vtk/MapsEllipsoidWithSource hippo_prototype.vtk ../../utilities/meshes/sphere320.vtk hippo1.vtk hippo2.vtk hippo3.vtk

# Creation of the prototype for the amygdala:
../../utilities/vtk/MapsEllipsoidWithSource amyg_prototype.vtk ../../utilities/meshes/sphere80.vtk amygdala1.vtk amygdala2.vtk amygdala3.vtk

