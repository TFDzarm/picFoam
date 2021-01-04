#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Set application name
application=$(getApplication)

#gmsh -3 Mesh1D.geo
#gmshToFoam Mesh1D.msh
#changeDictionary

blockMesh

picInitialise

runApplication $application
