#!/bin/bash
cp skull_australopithecus__t_19.vtk deformed_australopithecus.vtk
../../deformetrica/bin/ShootAndFlow2 paramDiffeos.xml -1 CP_final.txt Mom_final.txt paramCurves.xml deformed_australopithecus.vtk
../../deformetrica/bin/ShootAndFlow2 paramDiffeos.xml -1 CP_final.txt Mom_final.txt paramCurves.xml skull_neandertalis.vtk

