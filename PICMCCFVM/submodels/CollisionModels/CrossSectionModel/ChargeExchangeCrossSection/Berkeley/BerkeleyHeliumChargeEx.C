/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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

#include "BerkeleyHeliumChargeEx.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyHeliumChargeExCS<CloudType,Type>::BerkeleyHeliumChargeExCS
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
Foam::BerkeleyHeliumChargeExCS<CloudType,Type>::~BerkeleyHeliumChargeExCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyHeliumChargeExCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar crossx;

    if(eVEnergy < 0.01)//min
    {
        eVEnergy = 0.01;
    }
    if(eVEnergy < 377.8)//cutoff
    {
        crossx = 1.2996e-19  - 7.8872e-23*eVEnergy + 1.9873e-19/::sqrt(eVEnergy);
    }
    else
    {
        crossx = (1.5554e-18*::log(eVEnergy)/eVEnergy) + 8.6407e-20;
    }
    return crossx;
}

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyHeliumChargeExCS<CloudType,Type>::threshold() const
{
    return 0.0;
}

// ************************************************************************* //
