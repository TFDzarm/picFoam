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

#include "SentokuKempNonRelCorrection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SentokuKempNonRelCorrection<CloudType>::SentokuKempNonRelCorrection
(
    const dictionary& dict,
    CloudType& owner
)
:
    WeightCorrectionModel<CloudType>(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SentokuKempNonRelCorrection<CloudType>::~SentokuKempNonRelCorrection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::SentokuKempNonRelCorrection<CloudType>::correctVelocity(
                typename CloudType::parcelType* parcelP,
                typename CloudType::parcelType* parcelQ,
                vector preUp,
                vector preUq
                )
{
    if(!this->active())
        return;

    namespace cm = constant::mathematical;

    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    if(parcelP->nParticle() > parcelQ->nParticle())
    {
        // See SentokuKempCorrection.C for details
        scalar mass = cloud.constProps(parcelP->typeId()).mass();

        scalar wr = parcelQ->nParticle()/parcelP->nParticle();
        vector afterU = (1.0 - wr) * preUp + wr * parcelP->U();
        scalar afterE = (1.0 - wr) * 0.5 * mass * (preUp&preUp) + wr * 0.5 * mass * (parcelP->U()&parcelP->U());

        scalar deltaE = afterE-0.5*mass*(afterU&afterU);

        if(deltaE > 0.0)
        {
            vector e2(0.0,1.0,0.0);
            vector e3(0.0,0.0,1.0);
            scalar u_perp_mag = sqrt(sqr(afterU.y())+sqr(afterU.z()));
            if(u_perp_mag > VSMALL)
            {
                e2 = vector(0.0,afterU.z(),-afterU.y());
                e2 /= u_perp_mag;
                e3 = vector( sqr(u_perp_mag),
                             -(afterU.x()*afterU.y()),
                             -(afterU.x()*afterU.z()));
                e3 /= (mag(afterU) * u_perp_mag);
            }

            scalar eps = cm::twoPi*rndGen.scalar01();
            vector u_scat = afterU + sqrt(deltaE*2.0/mass) * (e2 * cos(eps) + e3 * sin(eps));
            parcelP->U() = u_scat;
        }
    }
    else if(parcelP->nParticle() < parcelQ->nParticle())
    {
        scalar mass = cloud.constProps(parcelQ->typeId()).mass();

        scalar wr = parcelP->nParticle()/parcelQ->nParticle();
        vector afterU = (1.0 - wr) *preUq + wr*parcelQ->U();
        scalar afterE = (1.0 - wr) *0.5*mass*(preUq&preUq) + wr*0.5*mass*(parcelQ->U()&parcelQ->U());
        scalar deltaE = afterE-0.5*mass*(afterU&afterU);

        if(deltaE > 0.0)
        {
            vector e2(0.0,1.0,0.0);
            vector e3(0.0,0.0,1.0);

            scalar u_perp_mag = sqrt(sqr(afterU.y())+sqr(afterU.z()));
            if(u_perp_mag > VSMALL)
            {
                e2 = vector(0.0,afterU.z(),-afterU.y());
                e2 /= u_perp_mag;
                e3 = vector( sqr(u_perp_mag),
                             -(afterU.x()*afterU.y()),
                             -(afterU.x()*afterU.z()));
                e3 /= (mag(afterU) * u_perp_mag);
            }
            scalar eps = cm::twoPi*rndGen.scalar01();
            vector u_scat = afterU + sqrt(deltaE*2.0/mass) * (e2 * cos(eps) + e3 * sin(eps));
            parcelQ->U() = u_scat;
        }
    }
}


// ************************************************************************* //
