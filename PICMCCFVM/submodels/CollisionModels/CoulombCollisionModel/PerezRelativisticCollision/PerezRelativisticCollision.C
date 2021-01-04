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

#include "PerezRelativisticCollision.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PerezRelativisticCollision<CloudType>::PerezRelativisticCollision
(
    const dictionary& dict,
    CloudType& cloud
)
:
    CoulombCollisionModel<CloudType>(dict,cloud,"CoulombCollision")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PerezRelativisticCollision<CloudType>::~PerezRelativisticCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PerezRelativisticCollision<CloudType>::active() const
{
    return true;
}

template<class CloudType>
void Foam::PerezRelativisticCollision<CloudType>::collide
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

    vector& UP = parcelP->U();
    vector& UQ = parcelQ->U();

    namespace cm = constant::mathematical;
    namespace ce = constant::electromagnetic;
    namespace cu = constant::universal;

    scalar deltaT = cloud.mesh().time().deltaTValue();

    //Get particle properties
    label idP = parcelP->typeId();
    scalar massP = cloud.constProps(idP).mass();
    scalar chargeP = parcelP->charge();

    label idQ = parcelQ->typeId();
    scalar massQ = cloud.constProps(idQ).mass();
    scalar chargeQ = parcelQ->charge();

    //Calculate the Lorentz factors
    scalar gammaP = 1.0/sqrt(1.0-sqr(mag(UP)/cu::c.value()));
    scalar gammaQ = 1.0/sqrt(1.0-sqr(mag(UQ)/cu::c.value()));

    //Relativistic momentum
    vector pP = massP*gammaP*UP;
    vector pQ = massQ*gammaQ*UQ;

    //Relative velocity
    vector vc_r = (pP + pQ) / (massP*gammaP+massQ*gammaQ);
    scalar magSqr_vcr = vc_r&vc_r;

    scalar gamma_c, gammac_SrqVcr, gammaP_CoM, gammaQ_CoM;
    vector pP_CoM/*,pQ_CoM*/;
    if(magSqr_vcr == 0.0) {
        gamma_c = 1.0;
        gammac_SrqVcr = 0.5;

        pP_CoM = pP;
        gammaP_CoM = gammaP;
        gammaQ_CoM = gammaQ;
    }
    else {
        //Transformation to the CoM frame of reference
        gamma_c = 1.0/sqrt( 1.0 - magSqr_vcr/sqr(cu::c.value()));
        gammac_SrqVcr = ((gamma_c -1.0)/magSqr_vcr);

        pP_CoM = pP + (gammac_SrqVcr*(vc_r&UP)-gamma_c)*massP*gammaP*vc_r;
        //pQ_CoM = pQ + (gammac_SrqVcr*(vc_r&parcelQ->U())-gamma_c)*massQ*gammaQ*vc_r;
        gammaP_CoM = (1.0-(vc_r&UP)/sqr(cu::c.value()))*gammaP*gamma_c;
        gammaQ_CoM = (1.0-(vc_r&UQ)/sqr(cu::c.value()))*gammaQ*gamma_c;
    }
    scalar mag_pP_CoM = mag(pP_CoM);
    if(mag_pP_CoM == 0.0)//Happens for some reason catch this...
        return;

    //Calculate or get the coulomb log
    scalar coulombLog;
    if(this->coulombLog_[idP][idQ] == 0.0) {
        //Impact parameter
        scalar b0 = chargeP*chargeQ/(4.0*cm::pi*ce::epsilon0.value()*sqr(cu::c.value()))*gamma_c/(massP*gammaP+massQ*gammaQ)*(1.0+massP*gammaP_CoM*massQ*gammaQ_CoM*sqr(cu::c.value())/sqr(mag_pP_CoM));

        scalar b_min = max(constant::universal::h.value()/(2*mag_pP_CoM),mag(b0));//choose minimum: de-Broglie wavelength or b0

        coulombLog = max(2.0,0.5*log(1.0+sqr(debyeLength)/sqr(b_min)));//FIXME: Performance !!!! this one only needs to calculated per cell -> move me !!!
    }
    else
        coulombLog = this->coulombLog_[idP][idQ];


    scalar s_term1 = gamma_c*mag_pP_CoM/(massP*gammaP+massQ*gammaQ);
    scalar s_term2 = sqr(1.0+massP*gammaP_CoM*massQ*gammaQ_CoM*sqr(cu::c.value())/sqr(mag_pP_CoM));
    scalar s = (nP*nQ/nPQ)*deltaT*coulombLog*sqr(chargeP*chargeQ)/(4.0*cm::pi*sqr(ce::epsilon0.value())*pow4(cu::c.value())*massP*gammaP*massQ*gammaQ)*s_term1*s_term2;

    //low temperatur correction
    scalar v_rel = mag_pP_CoM*(massP*gammaP+massQ*gammaQ)/(massP*gammaP_CoM*massQ*gammaQ_CoM*gamma_c);//should we use v1-v2 here??? Paper uses v1-v2 but talks about replacing with this, smilie uses this
    //scalar v_rel = mag(parcelP->U()-parcelQ->U());

    scalar s_max = pow(4.0*cm::pi/3.0,1.0/3.0) * (nP*nQ/nPQ)*deltaT*(massP+massQ)/max(massP*pow(nP,2.0/3.0),massQ*pow(nQ,2.0/3.0)) * v_rel;//max(mass1*::pow(n1,2.0/3.0),mass2*::pow(n2,2.0/3.0))
    s = min(s_max,s);

    //Scattering angles
    scalar cosChi = calculate_cosChi(s);
    scalar sinChi = sqrt( 1.0 - cosChi*cosChi );
    scalar eps = 2.0*cm::pi*cloud.rndGen().scalar01();

    scalar sinCcosE = sinChi*cos(eps);
    scalar sinCsinE = sinChi*sin(eps);

    //Perform the collision in the CoM frame
    vector new_pP_CoM;
    scalar p_perp = sqrt(pP_CoM.x()*pP_CoM.x() + pP_CoM.y()*pP_CoM.y());
    if( p_perp > 1.e-10*mag_pP_CoM ) {
        scalar inv_p_perp = 1.0/p_perp;
        new_pP_CoM = vector( (pP_CoM.x() * pP_CoM.z() * sinCcosE - pP_CoM.y() * mag_pP_CoM * sinCsinE) * inv_p_perp + pP_CoM.x() * cosChi   ,\
                             (pP_CoM.y() * pP_CoM.z() * sinCcosE + pP_CoM.x() * mag_pP_CoM * sinCsinE) * inv_p_perp + pP_CoM.y() * cosChi   ,\
                              -p_perp * sinCcosE  +  pP_CoM.z() * cosChi);
    }
    else
    {
        new_pP_CoM = vector(mag_pP_CoM*sinCcosE,mag_pP_CoM*sinCsinE,mag_pP_CoM*cosChi);
    }
    //Update the particles velocities
    vector new_pP = new_pP_CoM + vc_r*(gammac_SrqVcr*(vc_r&new_pP_CoM)+massP*gammaP_CoM*gamma_c);
    scalar new_gammaP = sqrt(1.0+sqr(mag(new_pP)/(massP*cu::c.value())));

    vector new_pQ = -1.0 * new_pP_CoM + vc_r*(gammac_SrqVcr*(vc_r&(-1.0*new_pP_CoM))+massQ*gammaQ_CoM*gamma_c);
    scalar new_gammaQ = sqrt(1.0+sqr(mag(new_pQ)/(massQ*cu::c.value())));

    vector preUP = parcelP->U();
    vector preUQ = parcelQ->U();

    UP = new_pP/(massP*new_gammaP);
    UQ = new_pQ/(massQ*new_gammaQ);

    this->weightCorrection().correctVelocity(parcelP,parcelQ,preUP,preUQ);//Weight correction
}

template<class CloudType>
scalar Foam::PerezRelativisticCollision<CloudType>::calculate_cosChi(scalar s)
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
