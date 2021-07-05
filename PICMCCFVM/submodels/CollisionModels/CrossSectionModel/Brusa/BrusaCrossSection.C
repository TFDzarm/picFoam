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

#include "BrusaCrossSection.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::BrusaCrossSectionBase::BrusaCrossSectionBase
(
    const dictionary& dict
)
:    s_(undefined)
{
    //Which dataset do we use?
    word species = dict.lookup("species");
    if(species == "Ne")
        s_ = Ne;
    else if(species == "Ar")
        s_ = Ar;
    else if(species == "Kr")
        s_ = Kr;
    else if(species == "Xe")
        s_ = Xe;
    else
        FatalErrorInFunction << "No elastic cross section for species " << species << " valid options are: (Ne, Ar, Kr, Xe)" << abort(FatalError);
}

template<class CloudType, Foam::crossSectionType Type>
Foam::BrusaCrossSection<CloudType,Type>::BrusaCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:   BrusaCrossSectionBase(dict.subDict(typeName + "Coeffs")),
    CrossSectionModel<CloudType,Type>(dict, cloud, typeName ,associatedTypeId)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BrusaCrossSectionBase::~BrusaCrossSectionBase()
{}


template<class CloudType, Foam::crossSectionType Type>
Foam::BrusaCrossSection<CloudType,Type>::~BrusaCrossSection()
{}


//Elastic Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::BrusaCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    BrusaCrossSectionBase(dict.subDict(typeName + "Coeffs")),
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronElasticCS>(dict, cloud, typeName ,associatedTypeId)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::~BrusaCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::crossSection(scalar eVEnergy) const
{
    //Cut of where the model is not valid
    if(eVEnergy < minEel_[s_])
        eVEnergy = minEel_[s_];

    if(eVEnergy > maxEel_[s_])
        eVEnergy=maxEel_[s_];//FIXME: Should this be done?

    eVEnergy/=1000.0;//in keV

    //Calculate the cross section according to Brusa
    scalar Qel = 1.0/(R1_[s_]*(S1_[s_]+eVEnergy)) + 1.0/(R2_[s_]*(S2_[s_]+eVEnergy));
    if(S2_[s_] > 0.0)
    {
        Qel += 2.0*::sqrt(S1_[s_]*S2_[s_]/(R1_[s_]*R2_[s_]))*::abs((::log((eVEnergy/S2_[s_]+1.0)/(eVEnergy/S1_[s_]+1.0))));
    }
    Qel*=1.e-20;
    return Qel;
}

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronElasticCS>::threshold() const
{
    return 0.0;
}



//Excitation Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::BrusaCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    BrusaCrossSectionBase(dict.subDict(typeName + "Coeffs")),
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronExciationCS>(dict, cloud, typeName ,associatedTypeId)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::~BrusaCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::crossSection(scalar eVEnergy) const
{
    //Cut of where the model is not valid
    if(eVEnergy < minEex_[s_])
        return 0.0;

    if(eVEnergy > maxEex_[s_])
        eVEnergy=maxEex_[s_];//FIXME: Should this be done?

    scalar EkeV = eVEnergy / 1000.0;

    //Calculate the cross section according to Brusa
    scalar Qex = 1.0/(F_[s_]*(G_[s_]+EkeV))*::log(eVEnergy/Eex_[s_]);
    Qex*=1.e-20;
    return Qex;
}

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronExciationCS>::threshold() const
{
    return minEex_[s_];
}


//Ionization Specialization
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::BrusaCrossSection
(
    const dictionary& dict,
    CloudType& cloud,
    const label& associatedTypeId
)
:
    BrusaCrossSectionBase(dict.subDict(typeName + "Coeffs")),
    CrossSectionModel<CloudType,Foam::crossSectionType::ElectronIonizationCS>(dict, cloud, typeName ,associatedTypeId)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::~BrusaCrossSection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::crossSection(scalar eVEnergy) const
{
    //Cut of where the model is not valid
    if(eVEnergy < minEion_[s_])
        return 0.0;

    if(eVEnergy > maxEion_[s_])
        eVEnergy=maxEion_[s_];//FIXME: Should this be done?


    scalar y = eVEnergy/Eion_[s_];

    eVEnergy/=1000.0;//in keV
    scalar x = eVEnergy/P_[s_];

    if(x < 1.0)
        x = 1.0;
    if(y < 1)
        y = 1.0;

    //Calculate the cross section according to Brusa
    scalar xySqrt = ::sqrt((y-1.0)/(x+1.0));
    scalar Qion = (L_[s_]/(M_[s_]+x)+N_[s_]/x)*xySqrt*xySqrt*xySqrt*(1.0+(2.0/3.0)*(1.0-1.0/(2.0*x))*::log(2.7+::sqrt(x-1.0)));
    Qion*=1.e-20;
    return Qion;
}

template<class CloudType>
Foam::scalar Foam::BrusaCrossSection<CloudType,Foam::crossSectionType::ElectronIonizationCS>::threshold() const
{
    return minEion_[s_];
}
// ************************************************************************* //
