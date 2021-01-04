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

#include "NanbuCollision.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuCollision<CloudType>::NanbuCollision
(
    const dictionary& dict,
    CloudType& cloud
)
:
    CoulombCollisionModel<CloudType>(dict,cloud,"CoulombCollision")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuCollision<CloudType>::~NanbuCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NanbuCollision<CloudType>::active() const
{
    return true;
}

template<class CloudType>
void Foam::NanbuCollision<CloudType>::collide
(
    typename CloudType::parcelType* parcelP,
    typename CloudType::parcelType* parcelQ,
    scalar debyeLength,
    scalar nP,
    scalar nQ,
    scalar nPQ
)
{
    CloudType& cloud(this->owner());

    namespace cm = constant::mathematical;
    namespace ce = constant::electromagnetic;
    namespace cu = constant::universal;

    //Get particle properties
    vector& UP = parcelP->U();
    vector& UQ = parcelQ->U();

    label typeIdP = parcelP->typeId();
    label typeIdQ = parcelQ->typeId();

    scalar deltaT = cloud.mesh().time().deltaTValue();

    scalar massP = cloud.constProps(typeIdP).mass();
    scalar chargeP = parcelP->charge();

    scalar massQ = cloud.constProps(typeIdQ).mass();
    scalar chargeQ = parcelQ->charge();

    scalar reMass = massP*massQ/(massP+massQ);

    vector g = UP-UQ;
    scalar g_mag = mag(g);

    //Impact parameter
    scalar b0 = chargeP*chargeQ/(2.0*cm::pi*ce::epsilon0.value()*reMass*g_mag*g_mag);//FIXME: Needs to use <g^2> instead if g^2...
    scalar b_min = max(cu::h.value()/(2.0*(g_mag*reMass)),b0);// Perez et al.: de-Broglie wavelength (unsure if g and reMass should be used), b0 (Nanbu97)


    //Get the coulomb log
    scalar coulombLog;
    if(this->coulombLog_[typeIdP][typeIdQ] == 0.0)
        coulombLog = max(2.0,0.5*::log(1.0+::pow(debyeLength,2)/::pow(b_min,2)));//FIXME: performance
    else
        coulombLog = this->coulombLog_[typeIdP][typeIdQ];

    //Calculate the scattering angles
    scalar s = coulombLog/(4.0*cm::pi) * ::pow(chargeP*chargeQ/(ce::epsilon0.value()*reMass),2)*(nQ*nP/nPQ)*::pow(g_mag,-3.0)*deltaT;
    scalar cosChi = calculate_cosChi(s);
    scalar sinChi = sqrt( 1.0 - cosChi*cosChi );

    scalar g_perp = ::sqrt(::pow(g.y(),2)+::pow(g.z(),2));
    scalar eps = 2.0*cm::pi*cloud.rndGen().scalar01();

    vector h(g_perp*::cos(eps),
             -1.0*(g.y()*g.x()*::cos(eps)+g_mag*g.z()*::sin(eps))/g_perp,
             -1.0*(g.z()*g.x()*::cos(eps)-g_mag*g.y()*::sin(eps))/g_perp);

    vector preUP = UP;
    vector preUQ = UQ;

    //Update the velocities
    UP = UP - massQ/(massP+massQ)*(g*(1.0-cosChi)+h*sinChi);
    UQ = UQ + massP/(massP+massQ)*(g*(1.0-cosChi)+h*sinChi);
    this->weightCorrection().correctVelocity(parcelP,parcelQ,preUP,preUQ);//Call weight correction
}

template<class CloudType>
scalar Foam::NanbuCollision<CloudType>::calculate_cosChi(scalar s)
{
    CloudType& cloud(this->owner());

    scalar A, invA;
    scalar U = cloud.rndGen().scalar01();

    if( s < 0.1 ) {
        if ( U<0.0001 ) U=0.0001; // ensures cos_chi > 0
        return 1.0 + s*::log(U);
    }
    if( s < 3.0  ) {
        // the polynomial has been modified from the article in order to have a better form // Snippet from Smilei
        invA = 0.00569578 +(0.95602 + (-0.508139 + (0.479139 + ( -0.12789 + 0.0238957*s )*s )*s )*s )*s;
        A = 1.0/invA;
        return  invA  * ::log( ::exp(-A) + 2.*U*::sinh(A) );
    }
    if( s < 6.0  ) {
        A = 3.0*exp(-s);
        return (1.0/A) * ::log( ::exp(-A) + 2.0*U*::sinh(A) );
    }
    return 2.0*U - 1.0;
}


// ************************************************************************* //
