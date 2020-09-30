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

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundaryEvent<CloudType>::BoundaryEvent(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    associatedPatches_(),
    interactionTable_(owner.mesh().boundaryMesh().size(),false)
{}

template<class CloudType>
Foam::BoundaryEvent<CloudType>::BoundaryEvent
(
    const dictionary& dict,
    CloudType& owner,
    const List<label>& associatedPatches
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dictionary::null),
    associatedPatches_(associatedPatches),
    interactionTable_(owner.mesh().boundaryMesh().size(),false)
{
    initInteractionTable();
}


template<class CloudType>
Foam::BoundaryEvent<CloudType>::BoundaryEvent
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const List<label>& associatedPatches
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    associatedPatches_(associatedPatches),
    interactionTable_(owner.mesh().boundaryMesh().size(),false)
{
    initInteractionTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundaryEvent<CloudType>::~BoundaryEvent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BoundaryEvent<CloudType>::initInteractionTable()
{
    forAll(associatedPatches_,i)
        interactionTable_[associatedPatches_[i]] = true;
}

template<class CloudType>
const CloudType& Foam::BoundaryEvent<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::BoundaryEvent<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::BoundaryEvent<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::BoundaryEvent<CloudType>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType>
const Foam::List<Foam::label>& Foam::BoundaryEvent<CloudType>::associatedPatches() const
{
    return associatedPatches_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoundaryEventModelNew.C"

// ************************************************************************* //
