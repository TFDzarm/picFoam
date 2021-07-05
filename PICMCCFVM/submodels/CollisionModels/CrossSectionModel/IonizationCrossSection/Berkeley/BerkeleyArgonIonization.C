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

#include "BerkeleyArgonIonization.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonIonizationCS<CloudType,Type>::BerkeleyArgonIonizationCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId),
    ionizationCS_()
{
    //Create interpolation table to speed up the calculation of the cross section
    label tableSize = 500;

    ionizationCS_.setSize(tableSize);

    for(label i = 0; i < tableSize; i++)
    {
        scalar eV = i;
        if(eV < threshold())
            ionizationCS_[i] = 0.0;
        else
            ionizationCS_[i] = (1.3596e-18/eV)*::log((eV +120.0/eV)/15.76)*(::atan((eV*eV -9.76*eV +2.4)/(20.6*eV +206)) + ::atan((2*eV -80.0)/(10.3*eV +103.0)));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonIonizationCS<CloudType,Type>::~BerkeleyArgonIonizationCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonIonizationCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qiz = 0.0;
    label i(eVEnergy);
    //Interpolate in table range
    if(i < ionizationCS_.size()-1)
    {
        scalar alpha = eVEnergy - i;
        Qiz = ionizationCS_[i]+alpha*(ionizationCS_[i+1]-ionizationCS_[i]);
    }
    else
        Qiz = ionizationCS_.last();
    return Qiz;
}


template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonIonizationCS<CloudType,Type>::threshold() const
{
    return 15.76;
}

// ************************************************************************* //
