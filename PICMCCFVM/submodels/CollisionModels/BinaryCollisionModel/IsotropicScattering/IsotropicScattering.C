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

#include "IsotropicScattering.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropicScattering<CloudType>::IsotropicScattering
(
    const dictionary& dict,
    CloudType& cloud
)
:
    BinaryCollisionModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropicScattering<CloudType>::~IsotropicScattering()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::IsotropicScattering<CloudType>::active() const
{
    return true;
}

template<class CloudType>
void Foam::IsotropicScattering<CloudType>::updateVelocity
(
    typename CloudType::parcelType& pP,
    typename CloudType::parcelType& pQ
)
{
    CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    //References to the particles velocities
    vector& UP = pP.U();
    vector& UQ = pQ.U();

    Random& rndGen(cloud.rndGen());

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cR = mag(UP - UQ);

    scalar cosTheta = 2.0*rndGen.scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();

    vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    //Update the velocities
    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
}

template<class CloudType>
void Foam::IsotropicScattering<CloudType>::updateVelocity
(
    typename CloudType::parcelType& pP,
    vector& Uq,
    label idQ
)
{
    CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();
    vector& UP = pP.U();

    Random& rndGen(cloud.rndGen());

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(idQ).mass();

    vector Ucm = (mP*UP + mQ*Uq)/(mP + mQ);

    scalar cR = mag(UP - Uq);

    scalar cosTheta = 2.0*rndGen.scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();

    vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );


    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
    Uq = Ucm - postCollisionRelU*mP/(mP + mQ);//Update the background velocity as well
}

// ************************************************************************* //
