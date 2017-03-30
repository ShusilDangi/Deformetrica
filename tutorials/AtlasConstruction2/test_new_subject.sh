#!/bin/bash

mkdir -p Results_atlas
mv *__t*.vtk Results_atlas
mv *final* Results_atlas

../../deformetrica/bin/sparseMatching3 paramDiffeosRegistration.xml paramHippo.xml hippo_prototype_template.vtk hippo4.vtk paramAmyg.xml amyg_prototype_template.vtk amygdala4.vtk

