/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "CrossSectionModel.H"

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionModel<CloudType,Type>::CrossSectionModel(CloudType& owner, const label& associatedTypeId)
:
  dict_(dictionary::null),
  owner_(owner),
  coeffDict_(dictionary::null),
  associatedTypeId_(associatedTypeId)
{}

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionModel<CloudType,Type>::CrossSectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const label& associatedTypeId
)
:
  dict_(dict),
  owner_(owner),
  coeffDict_(dictionary::null),
  associatedTypeId_(associatedTypeId)
{}

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionModel<CloudType,Type>::CrossSectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const label& associatedTypeId
)
:
  dict_(dict),
  owner_(owner),
  coeffDict_(dict.subDict(type + "Coeffs")),
  associatedTypeId_(associatedTypeId)
{}

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionModel<CloudType,Type>::~CrossSectionModel()
{}


template<class CloudType, Foam::crossSectionType Type>
const CloudType& Foam::CrossSectionModel<CloudType,Type>::owner() const
{
    return owner_;
}


template<class CloudType, Foam::crossSectionType Type>
CloudType& Foam::CrossSectionModel<CloudType,Type>::owner()
{
    return owner_;
}


template<class CloudType, Foam::crossSectionType Type>
const Foam::dictionary& Foam::CrossSectionModel<CloudType,Type>::dict() const
{
    return dict_;
}


template<class CloudType, Foam::crossSectionType Type>
const Foam::dictionary&
Foam::CrossSectionModel<CloudType,Type>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType, Foam::crossSectionType Type>
Foam::label
Foam::CrossSectionModel<CloudType,Type>::associatedTypeId() const
{
    return associatedTypeId_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CrossSectionModelNew.C"

// ************************************************************************* //
