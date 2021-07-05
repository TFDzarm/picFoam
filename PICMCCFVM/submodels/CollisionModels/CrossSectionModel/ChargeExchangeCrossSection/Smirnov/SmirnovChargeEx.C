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

#include "SmirnovChargeEx.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::SmirnovChargeExCS<CloudType, Type>::SmirnovChargeExCS
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType, Type>(dict, cloud, typeName, associatedTypeId),
    coeff1_(readScalar(this->coeffDict().lookup("coeff1"))),
    coeff2_(readScalar(this->coeffDict().lookup("coeff2")))
{
    Info << "        |= sigma(e) = 1e-20 * (" << coeff1_ << " - " << coeff2_ << " ln e)^2" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::SmirnovChargeExCS<CloudType,Type>::~SmirnovChargeExCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::SmirnovChargeExCS<CloudType,Type>::crossSection(scalar eVEnergy) const
{
    if(eVEnergy <= 0.0)
        return 0.0;

    scalar s = coeff1_-coeff2_*log(eVEnergy);
    return (s*s)*1e-20;
}


template<class CloudType, Foam::crossSectionType Type>
Foam::scalar Foam::SmirnovChargeExCS<CloudType,Type>::threshold() const
{
    return 0.0;
}

// ************************************************************************* //
