#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Set application name
application=$(getApplication)

gmsh -3 Mesh.geo
gmshToFoam Mesh.msh
changeDictionary

picInitialise

runApplication $application
