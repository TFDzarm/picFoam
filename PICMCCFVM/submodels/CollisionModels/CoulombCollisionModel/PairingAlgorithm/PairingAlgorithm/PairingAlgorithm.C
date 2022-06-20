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

#include "PairingAlgorithm.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairingAlgorithm<CloudType>::PairingAlgorithm(CloudType& cloud, CoulombCollisionModel<CloudType>& collisionModel)
:
    dict_(dictionary::null),
    cloud_(cloud),
    coulombCollisionModel_(collisionModel),
    coeffDict_(dictionary::null),
    nCells_(cloud.mesh().nCells())
{
    reduce(nCells_, sumOp<label>());
}


template<class CloudType>
Foam::PairingAlgorithm<CloudType>::PairingAlgorithm
(
    const dictionary& dict,
    CloudType& cloud,
    CoulombCollisionModel<CloudType>& collisionModel,
    const word& type
)
:
    dict_(dict),
    cloud_(cloud),
    coulombCollisionModel_(collisionModel),
    coeffDict_(dict.subOrEmptyDict(type + "Coeffs")),
    nCells_(cloud.mesh().nCells())
{
    reduce(nCells_, sumOp<label>());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairingAlgorithm<CloudType>::~PairingAlgorithm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::PairingAlgorithm<CloudType>::cloud() const
{
    return cloud_;
}


template<class CloudType>
CloudType&
Foam::PairingAlgorithm<CloudType>::cloud()
{
    return cloud_;
}

template<class CloudType>
const Foam::dictionary&
Foam::PairingAlgorithm<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::PairingAlgorithm<CloudType>::coeffDict() const
{
    return coeffDict_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PairingAlgorithmNew.C"

// ************************************************************************* //
