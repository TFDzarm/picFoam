/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "CenterInitialization.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CenterInitialization<CloudType>::CenterInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CenterInitialization<CloudType>::~CenterInitialization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CenterInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());

    const wordList species(this->coeffDict().lookup("species"));

    forAll(cloud.mesh().cells(), celli)
    {

        forAll(species, i)
        {
            const word& moleculeName(species[i]);
            label typeId(findIndex(cloud.typeIdList(), moleculeName));
            if(typeId == -1)
                FatalErrorInFunction << "Species " << moleculeName << " is not defined in cloud" << nl << FatalErrorInFunction;

            const typename CloudType::parcelType::constantProperties& cP =
            cloud.constProps(typeId);

            vector U = cloud.equipartitionLinearVelocity
            (
                this->temperatures_[typeId],
                cP.mass()
            );

            cloud.addNewParcel(cloud.mesh().C()[celli], celli, U, typeId);
        }
    }
}

// ************************************************************************* //

