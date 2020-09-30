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

#include "NanbuCorrection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuWeightCorrectionModel<CloudType>::NanbuWeightCorrectionModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    WeightCorrectionModel<CloudType>(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuWeightCorrectionModel<CloudType>::~NanbuWeightCorrectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::NanbuWeightCorrectionModel<CloudType>::correctVelocity(
                typename CloudType::parcelType* parcelP,
                typename CloudType::parcelType* parcelQ,
                vector preUp,
                vector preUq
                )
{
    if(!this->active())
        return;

    scalar rndU = this->owner().rndGen().scalar01();

    // with prop rndU <= p.nParticle()/q.nParticle() the collision occurs
    if(rndU > parcelP->nParticle()/parcelQ->nParticle())
        parcelQ->U() = preUq;

    // with prop rndU <= q.nParticle()/p.nParticle() the collisiom occurs
    if(rndU > parcelQ->nParticle()/parcelP->nParticle())
        parcelP->U() = preUp;
}


// ************************************************************************* //
