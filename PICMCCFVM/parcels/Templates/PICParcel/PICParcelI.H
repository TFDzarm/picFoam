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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::PICParcel<ParcelType>::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


template<class ParcelType>
inline Foam::PICParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    d_(readScalar(dict.lookup("diameter"))),
    omega_(readScalar(dict.lookup("omega"))),
    charge_(readScalar(dict.lookup("charge"))),
    solveMovement_(readBool(dict.lookup("solveMovement"))),
    nParticle_(readScalar(dict.lookup("nEquivalentParticles")))
{}


template<class ParcelType>
inline Foam::PICParcel<ParcelType>::PICParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const vector& U,
    const label typeId
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    U_(U),
    typeId_(typeId),
    chargeModifier_(1),
    charge_(0.0),
    nParticle_(1.0)
{}


template<class ParcelType>
inline Foam::PICParcel<ParcelType>::PICParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const vector& U,
    const label typeId
)
:
    ParcelType(mesh, position, celli),
    U_(U),
    typeId_(typeId),
    chargeModifier_(1),
    charge_(0.0),
    nParticle_(1.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::PICParcel<ParcelType>::constantProperties::mass() const
{
    return mass_;
}


template<class ParcelType>
inline Foam::scalar Foam::PICParcel<ParcelType>::constantProperties::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar
Foam::PICParcel<ParcelType>::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


template<class ParcelType>
inline Foam::scalar
Foam::PICParcel<ParcelType>::constantProperties::omega() const
{
    return omega_;
}

template<class ParcelType>
inline Foam::scalar
Foam::PICParcel<ParcelType>::constantProperties::charge() const
{
    return charge_;
}

template<class ParcelType>
inline bool
Foam::PICParcel<ParcelType>::constantProperties::solveMovement() const
{
    return solveMovement_;
}

template<class ParcelType>
inline scalar
Foam::PICParcel<ParcelType>::constantProperties::nParticle() const
{
    return nParticle_;
}

// * * * * * * * * * * PICParcel Member Functions  * * * * * * * * * * //

template<class ParcelType>
inline Foam::label Foam::PICParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline const Foam::vector& Foam::PICParcel<ParcelType>::U() const
{
    return U_;
}

template<class ParcelType>
inline Foam::vector& Foam::PICParcel<ParcelType>::U()
{
    return U_;
}

template<class ParcelType>
inline Foam::label Foam::PICParcel<ParcelType>::chargeModifier() const
{
    return chargeModifier_;
}

template<class ParcelType>
inline Foam::label Foam::PICParcel<ParcelType>::Zstar() const
{
    return chargeModifier_;
}

template<class ParcelType>
inline Foam::label& Foam::PICParcel<ParcelType>::chargeModifier()
{
    return chargeModifier_;
}

template<class ParcelType>
inline Foam::scalar Foam::PICParcel<ParcelType>::charge() const
{
    return charge_;
}

template<class ParcelType>
inline Foam::scalar& Foam::PICParcel<ParcelType>::charge()
{
    return charge_;
}

template<class ParcelType>
inline Foam::scalar Foam::PICParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}

template<class ParcelType>
inline Foam::scalar& Foam::PICParcel<ParcelType>::nParticle()
{
    return nParticle_;
}

template<class ParcelType>
inline vector
Foam::PICParcel<ParcelType>::trackingData::boundaryVelocity() const
{
    return boundaryVelocity_;
}

template<class ParcelType>
inline vector&
Foam::PICParcel<ParcelType>::trackingData::boundaryVelocity()
{
    return boundaryVelocity_;
}

template<class ParcelType>
inline bool
Foam::PICParcel<ParcelType>::trackingData::requireResync() const
{
    return requireResync_;
}


template<class ParcelType>
inline bool&
Foam::PICParcel<ParcelType>::trackingData::requireResync()
{
    return requireResync_;
}

// ************************************************************************* //
