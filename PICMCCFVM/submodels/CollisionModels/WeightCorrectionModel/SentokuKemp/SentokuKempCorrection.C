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

#include "SentokuKempCorrection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SentokuKempCorrection<CloudType>::SentokuKempCorrection
(
    const dictionary& dict,
    CloudType& owner
)
:
    WeightCorrectionModel<CloudType>(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SentokuKempCorrection<CloudType>::~SentokuKempCorrection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::SentokuKempCorrection<CloudType>::correctVelocity(
                typename CloudType::parcelType* parcelP,
                typename CloudType::parcelType* parcelQ,
                vector preUp,
                vector preUq
                )
{
    if(!this->active())
        return;

    namespace cm = constant::mathematical;
    namespace cu = constant::universal;

    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    if(parcelP->nParticle() > parcelQ->nParticle())//first case P has more particles than Q
    {
        scalar mass = cloud.constProps(parcelP->typeId()).mass();
        scalar deltaW = parcelQ->nParticle()/parcelP->nParticle();

        //Lorentz factors
        scalar pre_gamma = 1.0/sqrt(1.0-sqr(mag(preUp)/cu::c.value()));
        scalar post_gamma = 1.0/sqrt(1.0-sqr(mag(parcelP->U())/cu::c.value()));

        //Momentum
        vector pre_pP = mass*pre_gamma*preUp;
        vector post_pP = mass*post_gamma*parcelP->U();

        //Energy
        scalar pre_eP = cu::c.value()*sqrt((pre_pP&pre_pP) + sqr(mass * cu::c.value()));//Total Energy
        scalar post_eP = cu::c.value()*sqrt((post_pP&post_pP) + sqr(mass * cu::c.value()));

        //Corrected momentum and energy
        scalar afterColE = (1.0 - deltaW) * pre_eP + deltaW * post_eP;
        vector afterColMom = (1.0 - deltaW) * pre_pP + deltaW * post_pP;

        scalar afterColMom_mag = mag(afterColMom);

        //Not 1+ afterColE ...! Because afterColE is the total energy
        scalar gamma_e = afterColE / (mass * sqr(cu::c.value()));
        scalar gamma_p = sqrt(1.0 + sqr(afterColMom_mag/(mass*cu::c.value())));

        //This should be the case, we correct this...
        if(gamma_p < gamma_e)
        {
            scalar delta_p = mass * cu::c.value() * sqrt(sqr(gamma_e)-sqr(gamma_p));
            scalar p_perp_mag = sqrt(sqr(afterColMom.y())+sqr(afterColMom.z()));

            //new coords as in EPOCH
            //default values
            vector e2(0.0,1.0,0.0);
            vector e3(0.0,0.0,1.0);

            if(p_perp_mag > VSMALL)//if vector is x direction skip/use defaults
            {
                e2 = vector(0.0,afterColMom.z(),-afterColMom.y());//first orthogonal vector to afterColMom
                e2 = e2 / p_perp_mag;
                e3 = vector( sqr(p_perp_mag),
                             -(afterColMom.x()*afterColMom.y()),
                             -(afterColMom.x()*afterColMom.z()));//second orthogonal found by cross product
                e3 /= (afterColMom_mag * p_perp_mag);
            }
            scalar eps = cm::twoPi*rndGen.scalar01();//choose random vector in orthogonal plane

            vector p_scat = afterColMom + delta_p * (e2 * cos(eps) + e3 * sin(eps));//random orthogonal direction
            scalar gamma_scat = sqrt(1.0+sqr(mag(p_scat)/(mass*cu::c.value())));
            parcelP->U() = p_scat/(mass*gamma_scat);//apply correction
        }
    }
    else if(parcelP->nParticle() < parcelQ->nParticle())//first case Q has more particles than P
    {
        scalar mass = cloud.constProps(parcelQ->typeId()).mass();

        scalar deltaW = parcelP->nParticle()/parcelQ->nParticle();

        scalar pre_gamma = 1.0/sqrt(1.0-sqr(mag(preUq)/cu::c.value()));
        scalar post_gamma = 1.0/sqrt(1.0-sqr(mag(parcelQ->U())/cu::c.value()));

        vector pre_pQ = mass*pre_gamma*preUq;
        vector post_pQ = mass*post_gamma*parcelQ->U();

        scalar pre_eQ = cu::c.value()*sqrt((pre_pQ&pre_pQ) + sqr(mass * cu::c.value()));//Total Energy
        scalar post_eQ = cu::c.value()*sqrt((post_pQ&post_pQ) + sqr(mass * cu::c.value()));

        scalar afterColE = (1.0 - deltaW) * pre_eQ + deltaW * post_eQ;
        vector afterColMom = (1.0 - deltaW) * pre_pQ + deltaW * post_pQ;

        scalar afterColMom_mag = mag(afterColMom);

        scalar gamma_e = afterColE / (mass * sqr(cu::c.value()));
        scalar gamma_p = sqrt(1.0 + sqr(afterColMom_mag/(mass*cu::c.value())));
        if(gamma_p < gamma_e)
        {
            scalar delta_p = mass * cu::c.value() * sqrt(sqr(gamma_e)-sqr(gamma_p));
            scalar p_perp_mag = sqrt(sqr(afterColMom.y())+sqr(afterColMom.z()));

            //new coords as in EPOCH
            //default values
            vector e2(0.0,1.0,0.0);
            vector e3(0.0,0.0,1.0);

            if(p_perp_mag > VSMALL) //first orthogonal vector to afterColMom
            {
                e2 = vector(0.0,afterColMom.z(),-afterColMom.y());//first orthogonal vector to afterColMom
                e2 = e2 / p_perp_mag;
                e3 = vector( sqr(p_perp_mag),
                             -(afterColMom.x()*afterColMom.y()),
                             -(afterColMom.x()*afterColMom.z()));//second orthogonal found by cross product
                e3 /= (afterColMom_mag * p_perp_mag);
            }
            scalar eps = cm::twoPi*rndGen.scalar01();//choose random vector in orthogonal plane

            vector p_scat = afterColMom + delta_p * (e2 * cos(eps) + e3 * sin(eps));//random orthogonal direction
            scalar gamma_scat = sqrt(1.0+sqr(mag(p_scat)/(mass*cu::c.value())));
            parcelQ->U() = p_scat/(mass*gamma_scat);//apply correction
        }
    }

}


// ************************************************************************* //
