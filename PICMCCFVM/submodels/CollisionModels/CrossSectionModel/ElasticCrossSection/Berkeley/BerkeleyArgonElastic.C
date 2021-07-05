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

#include "BerkeleyArgonElastic.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonElasticCS<CloudType,Type>::BerkeleyArgonElasticCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Type>(cloud, associatedTypeId),
    elasticLowEnergyCS_(),
    elasticCS_()
{
    //Create interpolation table to speed up the calculation of the cross section
    label tableSizeLowEnergy = 101;
    label tableSize = 500;

    elasticLowEnergyCS_.setSize(tableSizeLowEnergy);
    elasticCS_.setSize(tableSize);

    for(label i = 0; i < tableSizeLowEnergy; i++)
    {
        scalar eV = 0.01*i;
        if(eV < 0.2)
            elasticLowEnergyCS_[i] = 1.0/::pow(10.0,19.0 + eV/0.11);
        else
            elasticLowEnergyCS_[i] = 9.07e-19*::pow(eV,1.55)*::pow(eV+70.0,1.1)/::pow(14.0+eV,3.25);
    }
    for(label i = 0; i < tableSize; i++)
    {
        scalar eV = i;
        elasticCS_[i] = 9.07e-19*::pow(eV, 1.55)*::pow(eV+70.0, 1.10)/::pow(14.+eV, 3.25);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::BerkeleyArgonElasticCS<CloudType,Type>::~BerkeleyArgonElasticCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonElasticCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    scalar Qel = 0.0;
    //Interpolation for low energies
    if(eVEnergy < 1.0)
    {
        label i(100*eVEnergy);
        scalar alpha = 100*eVEnergy - i;
        Qel = elasticLowEnergyCS_[i] + alpha*(elasticLowEnergyCS_[i+1]-elasticLowEnergyCS_[i]);
    }
    else//Interpolation for high energies
    {
        label i(eVEnergy);
        if(i < elasticCS_.size()-1)
        {
            scalar alpha = eVEnergy-i;
            Qel = elasticCS_[i]+alpha*(elasticCS_[i+1]-elasticCS_[i]);
        }
        else
            Qel = elasticCS_.last();
    }
    return Qel;
}


template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::BerkeleyArgonElasticCS<CloudType,Type>::threshold() const
{
    return 0.0;
}

// ************************************************************************* //
