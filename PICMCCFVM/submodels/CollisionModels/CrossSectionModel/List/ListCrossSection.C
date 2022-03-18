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

#include "ListCrossSection.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

//Helper function, was removed in OF-9
template<class Type>
Type interpolateXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    label n = xOld.size();

    label lo = 0;
    for (lo=0; lo<n && xOld[lo]>x; ++lo)
    {}

    label low = lo;
    if (low < n)
    {
        for (label i=low; i<n; ++i)
        {
            if (xOld[i] > xOld[lo] && xOld[i] <= x)
            {
                lo = i;
            }
        }
    }

    label hi = 0;
    for (hi=0; hi<n && xOld[hi]<x; ++hi)
    {}

    label high = hi;
    if (high < n)
    {
        for (label i=high; i<n; ++i)
        {
            if (xOld[i] < xOld[hi] && xOld[i] >= x)
            {
                hi = i;
            }
        }
    }


    if (lo<n && hi<n && lo != hi)
    {
        return yOld[lo]
            + ((x - xOld[lo])/(xOld[hi] - xOld[lo]))*(yOld[hi] - yOld[lo]);
    }
    else if (lo == hi)
    {
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi];
    }
    else
    {
        return yOld[lo];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::ListCrossSection<CloudType,Type>::ListCrossSection
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
Foam::ListCrossSection<CloudType,Type>::~ListCrossSection()
{}


//Elastic Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::ListCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronElasticCS>(dict, cloud, "Elastic" + typeName ,associatedTypeId),
    energies_(this->coeffDict().lookup("energies")),
    values_(this->coeffDict().lookup("values"))
{
      if(energies_.size() != values_.size())
        FatalErrorInFunction << "Lists \"energies\" (size: " << energies_.size() << ") and \"values\" (size: " << values_.size() << ") are not of the same size" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::~ListCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::crossSection(scalar eVEnergy) const
{
    return interpolateXY(eVEnergy,energies_,values_);
}

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::threshold() const
{
    return 0.0;
}



//Excitation Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::ListCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronExciationCS>(dict, cloud, "Excitation" + typeName ,associatedTypeId),
    energies_(this->coeffDict().lookup("energies")),
    values_(this->coeffDict().lookup("values")),
    threshold_(readScalar(this->coeffDict().lookup("threshold")))
{
      if(energies_.size() != values_.size())
        FatalErrorInFunction << "Lists \"energies\" (size: " << energies_.size() << ") and \"values\" (size: " << values_.size() << ") are not of the same size" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::~ListCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::crossSection(scalar eVEnergy) const
{
    return interpolateXY(eVEnergy,energies_,values_);
}

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::threshold() const
{
    return threshold_;
}


//Ionization Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::ListCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronIonizationCS>(dict, cloud, "Ionization" + typeName ,associatedTypeId),
    energies_(this->coeffDict().lookup("energies")),
    values_(this->coeffDict().lookup("values")),
    threshold_(readScalar(this->coeffDict().lookup("threshold")))
{
      if(energies_.size() != values_.size())
        FatalErrorInFunction << "Lists \"energies\" (size: " << energies_.size() << ") and \"values\" (size: " << values_.size() << ") are not of the same size" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::~ListCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::crossSection(scalar eVEnergy) const
{
    return interpolateXY(eVEnergy,energies_,values_);
}

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::threshold() const
{
    return threshold_;
}



//Ion Elastic Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonElasticCS>::ListCrossSection
    (
        const dictionary& dict,
        CloudType& cloud,
        const label& associatedTypeId
        )
    :
    CrossSectionModel<CloudType,Foam::crossSectionType::IonElasticCS>(dict, cloud, "Elastic" + typeName ,associatedTypeId),
    energies_(this->coeffDict().lookup("energies")),
    values_(this->coeffDict().lookup("values"))
{
      if(energies_.size() != values_.size())
        FatalErrorInFunction << "Lists \"energies\" (size: " << energies_.size() << ") and \"values\" (size: " << values_.size() << ") are not of the same size" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonElasticCS>::~ListCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonElasticCS>::crossSection(scalar eVEnergy) const
{
    return interpolateXY(eVEnergy,energies_,values_);
}

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonElasticCS>::threshold() const
{
    return 0.0;
}




//Ion ChargeEx Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonChargeExCS>::ListCrossSection
    (
        const dictionary& dict,
        CloudType& cloud,
        const label& associatedTypeId
        )
    :
    CrossSectionModel<CloudType,Foam::crossSectionType::IonChargeExCS>(dict, cloud, "ChargeExchange" + typeName ,associatedTypeId),
    energies_(this->coeffDict().lookup("energies")),
    values_(this->coeffDict().lookup("value")),
    threshold_(readScalar(this->coeffDict().lookup("threshold")))
{
      if(energies_.size() != values_.size())
        FatalErrorInFunction << "Lists \"energies\" (size: " << energies_.size() << ") and \"values\" (size: " << values_.size() << ") are not of the same size" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonChargeExCS>::~ListCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonChargeExCS>::crossSection(scalar eVEnergy) const
{
    return interpolateXY(eVEnergy,energies_,values_);
}

template<class CloudType>
Foam::scalar Foam::ListCrossSection<CloudType,Foam::crossSectionType::IonChargeExCS>::threshold() const
{
    return threshold_;
}


// ************************************************************************* //
