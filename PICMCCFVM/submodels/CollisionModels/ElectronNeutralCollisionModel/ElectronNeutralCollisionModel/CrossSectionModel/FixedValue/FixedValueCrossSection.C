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

#include "FixedValueCrossSection.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::FixedValueCrossSection<CloudType,Type>::FixedValueCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:   CrossSectionModel<CloudType,Type>(dict, cloud, typeName ,associatedTypeId)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::FixedValueCrossSection<CloudType,Type>::~FixedValueCrossSection()
{}


//Elastic Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::FixedValueCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronElasticCS>(dict, cloud, "Elastic" + typeName ,associatedTypeId),
    value_(readScalar(this->coeffDict().lookup("value")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::~FixedValueCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::crossSection(scalar eVEnergy) const
{
    return value_;
}

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::threshold() const
{
    return 0.0;
}



//Excitation Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::FixedValueCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronExciationCS>(dict, cloud, "Exciation" + typeName ,associatedTypeId),
    value_(readScalar(this->coeffDict().lookup("value"))),
    threshold_(readScalar(this->coeffDict().lookup("value")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::~FixedValueCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::crossSection(scalar eVEnergy) const
{
    return value_;
}

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::threshold() const
{
    return threshold_;
}


//Ionization Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::FixedValueCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronIonizationCS>(dict, cloud, "Ionization" + typeName ,associatedTypeId),
    value_(readScalar(this->coeffDict().lookup("value"))),
    threshold_(readScalar(this->coeffDict().lookup("value")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::~FixedValueCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::crossSection(scalar eVEnergy) const
{
    return value_;//return the fixed value
}

template<class CloudType>
Foam::scalar Foam::FixedValueCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::threshold() const
{
    return threshold_;
}
// ************************************************************************* //
