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

#include "CrossSectionInitialization.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CrossSectionInitialization<CloudType>::CrossSectionInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName)
{
    bool iCM = readBool(dict.lookup("initalizeCollisionModels"));
    bool sME = readBool(dict.lookup("solveMaxwellEquations"));
    bool iLF = readBool(dict.lookup("initalizeLeapFrog"));
    if(!iCM) {
        FatalErrorInFunction << "initalizeCollisionModels has to be set to true" << abort(FatalError);
    }
    if(!sME) {
        FatalErrorInFunction << "solveMaxwellEquations has to be set to false" << abort(FatalError);
    }
    if(!iLF) {
        FatalErrorInFunction << "initalizeLeapFrog has to be set to false" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CrossSectionInitialization<CloudType>::~CrossSectionInitialization()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::CrossSectionInitialization<CloudType>::initialiseParticles()
{}
// ************************************************************************* //

