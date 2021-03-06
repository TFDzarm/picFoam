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

#include "BerkeleyHeliumElastic.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// ------------------------------- Electron-Neutral -------------------------------

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyHeliumElasticCS<CloudType,Type>::BerkeleyHeliumElasticCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyHeliumElasticCS<CloudType,Type>::~BerkeleyHeliumElasticCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyHeliumElasticCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qel = 8.5e-19/(pow(eVEnergy+10.0, 1.1));
    return Qel;
}


template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyHeliumElasticCS<CloudType,Type>::threshold() const
{
    return 0.0;
}

// ------------------------------- Ion-Neutral -------------------------------

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BerkeleyHeliumElasticCS<CloudType,crossSectionType::IonElasticCS>::BerkeleyHeliumElasticCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,crossSectionType::IonElasticCS>(cloud, associatedTypeId)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BerkeleyHeliumElasticCS<CloudType,crossSectionType::IonElasticCS>::~BerkeleyHeliumElasticCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BerkeleyHeliumElasticCS<CloudType,crossSectionType::IonElasticCS>::crossSection(scalar eVEnergy) const
{
    scalar crossx;
    if(eVEnergy < 0.01)
    {
        eVEnergy = 0.01;
    }

    crossx = 3.6463e-19/::sqrt(eVEnergy) - 7.9897e-21;
    if(crossx < 0)
    {
        return 0.0;
    }
    return crossx;
}

template<class CloudType>
Foam::scalar Foam::BerkeleyHeliumElasticCS<CloudType,crossSectionType::IonElasticCS>::threshold() const
{
    return 0.0;
}

// ************************************************************************* //
