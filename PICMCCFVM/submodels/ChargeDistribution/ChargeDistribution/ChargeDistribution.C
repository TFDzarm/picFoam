/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "ChargeDistribution.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeDistribution<CloudType>::ChargeDistribution(const word& fieldName, CloudType& owner)
:
    dict_(dictionary::null),
    cloud_(owner),
    coeffDict_(dictionary::null),
    field_(
        IOobject
        (
            fieldName,
            cloud_.mesh().time().timeName(),
            cloud_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        cloud_.mesh()    
    )
{}


template<class CloudType>
Foam::ChargeDistribution<CloudType>::ChargeDistribution
(
    const word& fieldName,
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    cloud_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    field_(
        IOobject
        (
            fieldName,
            cloud_.mesh().time().timeName(),
            cloud_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        cloud_.mesh()    
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeDistribution<CloudType>::~ChargeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::ChargeDistribution<CloudType>::cloud() const
{
    return cloud_;
}


template<class CloudType>
CloudType& Foam::ChargeDistribution<CloudType>::cloud()
{
    return cloud_;
}


template<class CloudType>
const Foam::dictionary& Foam::ChargeDistribution<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::ChargeDistribution<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
void Foam::ChargeDistribution<CloudType>::update()
{}

template<class CloudType>
void Foam::ChargeDistribution<CloudType>::reset()
{
    field_ = dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0,1,0), 0.0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChargeDistributionNew.C"

// ************************************************************************* //
