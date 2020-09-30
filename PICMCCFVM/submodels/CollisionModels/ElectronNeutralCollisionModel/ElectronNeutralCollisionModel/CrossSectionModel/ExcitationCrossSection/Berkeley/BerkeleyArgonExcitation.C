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

#include "BerkeleyArgonExcitation.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonExcitationCS<CloudType,Type>::BerkeleyArgonExcitationCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId),
    excitationCS_()
{
    label tableSize = 500;

    excitationCS_.setSize(tableSize);

    for(label i = 0; i < tableSize; i++)
    {
        scalar eV = i;
        if(eV < threshold())
            excitationCS_[i] = 0.0;
        else
            excitationCS_[i] = (3.85116e-19*::log(eV/3.4015) -4.85227e-19)/eV;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonExcitationCS<CloudType,Type>::~BerkeleyArgonExcitationCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonExcitationCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qex = 0.0;
    label i(eVEnergy);
    if(i < excitationCS_.size()-1)
    {
        scalar alpha = eVEnergy - i;
        Qex = excitationCS_[i]+alpha*(excitationCS_[i+1]-excitationCS_[i]);
    }
    else
        Qex = excitationCS_.last();
    return Qex;
}


template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonExcitationCS<CloudType,Type>::threshold() const
{
    return 11.55;
}

// ************************************************************************* //
