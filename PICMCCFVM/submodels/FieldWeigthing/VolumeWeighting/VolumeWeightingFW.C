/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 picFoam
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "VolumeWeightingFW.H"
#include "polyMeshTetDecomposition.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VolumeWeightingFW::VolumeWeightingFW
(
    const word& fieldName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    FieldWeigthing(fieldName, dict, mesh),
    Einterp_(autoPtr<interpolationCellPoint<vector>>(new interpolationCellPoint<vector>(field_)))
{}


Foam::VolumeWeightingFW::VolumeWeightingFW
(
    const VolumeWeightingFW& am
)
:
    FieldWeigthing(am),
    Einterp_(am.Einterp_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VolumeWeightingFW::~VolumeWeightingFW()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::VolumeWeightingFW::getFieldVector
(
    const barycentric& coordinates,
    const tetIndices& tetIs
) const
{
    return Einterp_->interpolate(coordinates,tetIs);
}


void Foam::VolumeWeightingFW::update()
{
    //This interpolates the field to the verticies and communicates with all processors
    //FIXME: Optimize? This calculates weights new every timestep. We do not support moving meshes...
    Einterp_.reset(new interpolationCellPoint<vector>(field_));
}

// ************************************************************************* //
