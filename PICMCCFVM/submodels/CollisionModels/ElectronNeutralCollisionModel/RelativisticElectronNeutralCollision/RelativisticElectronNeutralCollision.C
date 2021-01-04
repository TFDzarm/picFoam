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

#include "RelativisticElectronNeutralCollision.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RelativisticElectronNeutralCollision<CloudType>::RelativisticElectronNeutralCollision
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ElectronNeutralCollisionModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RelativisticElectronNeutralCollision<CloudType>::~RelativisticElectronNeutralCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::RelativisticElectronNeutralCollision<CloudType>::active() const
{
    return true;
}

template<class CloudType>
void Foam::RelativisticElectronNeutralCollision<CloudType>::updateVelocity(scalar eV, typename CloudType::parcelType& parcelE, typename CloudType::parcelType& parcelN)
{
    namespace cm = constant::mathematical;
    namespace cu = constant::universal;

    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    //Velocity references
    vector& Ue = parcelE.U();
    vector& Un = parcelN.U();

    //TypeId and mass
    label idE = parcelE.typeId();
    scalar massE = cloud.constProps(idE).mass();

    label idN = parcelN.typeId();
    scalar massN = cloud.constProps(idN).mass();

    //Lorentz factors for parcel E and N
    scalar gammaE = 1.0/sqrt(1.0-sqr(mag(Ue)/cu::c.value()));
    scalar gammaN = 1.0/sqrt(1.0-sqr(mag(Un)/cu::c.value()));

    //Relativstic momentum
    vector pE = massE*gammaE*Ue;
    vector pN = massN*gammaN*Un;

    //Relative veloctiy
    vector vc_r = (pE + pN) / (massE*gammaE+massN*gammaN);
    scalar magSqr_vcr = vc_r&vc_r;

    scalar gamma_c, gammac_SrqVcr, gammaE_CoM, gammaN_CoM;
    vector pE_CoM/*,pQ_CoM*/;
    if(magSqr_vcr == 0.0) {
        gamma_c = 1.0;
        gammac_SrqVcr = 0.5;

        pE_CoM = pE;
        gammaE_CoM = gammaE;
        gammaN_CoM = gammaN;
    }
    else {
        gamma_c = 1.0/sqrt( 1.0 - magSqr_vcr/sqr(cu::c.value()));
        gammac_SrqVcr = ((gamma_c -1.0)/magSqr_vcr);

        //Lorentz transformation to the Center of Mass system
        pE_CoM = pE + (gammac_SrqVcr*(vc_r&Ue)-gamma_c)*massE*gammaE*vc_r;
        //pQ_CoM = pQ + (gammac_SrqVcr*(vc_r&parcelQ->U())-gamma_c)*massQ*gammaQ*vc_r;
        gammaE_CoM = (1.0-(vc_r&Ue)/sqr(cu::c.value()))*gammaE*gamma_c;
        gammaN_CoM = (1.0-(vc_r&Un)/sqr(cu::c.value()))*gammaN*gamma_c;
    }
    scalar mag_pE_CoM = mag(pE_CoM);
    if(mag_pE_CoM == 0.0)//? rarely happens
        return;

    scalar cosChi, sinChi;
    if(eV < 1e-30) {
        cosChi = 1.0;
        sinChi = 0.0;
    }
    else {
        cosChi = 1.0/eV*(eV+2.-2.*::pow(1.+eV,rndGen.scalar01()));
        if(cosChi*cosChi > 1.0) {//floating point errors this happens rarely...
            cosChi = 1.0;
            sinChi = 0.0;
        }
        else
            sinChi = sqrt( 1.0 - cosChi*cosChi );
    }

    scalar eps = cm::twoPi*rndGen.scalar01();

    scalar sinCcosE = sinChi*cos(eps);
    scalar sinCsinE = sinChi*sin(eps);

    vector new_pE_CoM;
    scalar p_perp = sqrt(pE_CoM.x()*pE_CoM.x() + pE_CoM.y()*pE_CoM.y());
    if( p_perp > 1.e-10*mag_pE_CoM ) {
        //Scattering in the CoM frame
        scalar inv_p_perp = 1.0/p_perp;
        new_pE_CoM = vector( (pE_CoM.x() * pE_CoM.z() * sinCcosE - pE_CoM.y() * mag_pE_CoM * sinCsinE) * inv_p_perp + pE_CoM.x() * cosChi   ,\
                             (pE_CoM.y() * pE_CoM.z() * sinCcosE + pE_CoM.x() * mag_pE_CoM * sinCsinE) * inv_p_perp + pE_CoM.y() * cosChi   ,\
                              -p_perp * sinCcosE  +  pE_CoM.z() * cosChi);
    }
    else {
        new_pE_CoM = vector(mag_pE_CoM*sinCcosE,mag_pE_CoM*sinCsinE,mag_pE_CoM*cosChi);
    }

    //Transformation back to the parcel's frame
    vector new_pE = new_pE_CoM + vc_r*(gammac_SrqVcr*(vc_r&new_pE_CoM)+massE*gammaE_CoM*gamma_c);
    scalar new_gammaE = sqrt(1.0+sqr(mag(new_pE)/(massE*cu::c.value())));

    vector new_pN = -1.0 * new_pE_CoM + vc_r*(gammac_SrqVcr*(vc_r&(-1.0*new_pE_CoM))+massN*gammaN_CoM*gamma_c);
    scalar new_gammaN = sqrt(1.0+sqr(mag(new_pN)/(massN*cu::c.value())));

    Ue = new_pE/(massE*new_gammaE);
    Un = new_pN/(massN*new_gammaN);
}

//Same as above, but uses sampled velocity Un to calculate the scattering
template<class CloudType>
void Foam::RelativisticElectronNeutralCollision<CloudType>::updateVelocity(scalar eV, typename CloudType::parcelType& parcelE, vector& Un, label idN)
{
    namespace cm = constant::mathematical;
    namespace cu = constant::universal;

    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    vector& Ue = parcelE.U();

    label idE = parcelE.typeId();
    scalar massE = cloud.constProps(idE).mass();

    scalar massN = cloud.constProps(idN).mass();

    scalar gammaE = 1.0/sqrt(1.0-sqr(mag(Ue)/cu::c.value()));
    scalar gammaN = 1.0/sqrt(1.0-sqr(mag(Un)/cu::c.value()));

    vector pE = massE*gammaE*Ue;
    vector pN = massN*gammaN*Un;

    vector vc_r = (pE + pN) / (massE*gammaE+massN*gammaN);
    scalar magSqr_vcr = vc_r&vc_r;

    scalar gamma_c, gammac_SrqVcr, gammaE_CoM, gammaN_CoM;
    vector pE_CoM/*,pQ_CoM*/;
    if(magSqr_vcr == 0.0) {
        gamma_c = 1.0;
        gammac_SrqVcr = 0.5;

        pE_CoM = pE;
        gammaE_CoM = gammaE;
        gammaN_CoM = gammaN;
    }
    else {
        gamma_c = 1.0/sqrt( 1.0 - magSqr_vcr/sqr(cu::c.value()));
        gammac_SrqVcr = ((gamma_c -1.0)/magSqr_vcr);

        pE_CoM = pE + (gammac_SrqVcr*(vc_r&Ue)-gamma_c)*massE*gammaE*vc_r;
        //pQ_CoM = pQ + (gammac_SrqVcr*(vc_r&parcelQ->U())-gamma_c)*massQ*gammaQ*vc_r;
        gammaE_CoM = (1.0-(vc_r&Ue)/sqr(cu::c.value()))*gammaE*gamma_c;
        gammaN_CoM = (1.0-(vc_r&Un)/sqr(cu::c.value()))*gammaN*gamma_c;
    }
    scalar mag_pE_CoM = mag(pE_CoM);
    if(mag_pE_CoM == 0.0)//? rarely happens
        return;

    scalar cosChi;
    if(eV < 1e-30)
        cosChi = 1.0;
    else
        cosChi = 1.0/eV*(eV+2.-2.*::pow(1.+eV,rndGen.scalar01()));

    scalar sinChi = sqrt( 1.0 - cosChi*cosChi );
    scalar eps = cm::twoPi*rndGen.scalar01();

    scalar sinCcosE = sinChi*cos(eps);
    scalar sinCsinE = sinChi*sin(eps);

    vector new_pE_CoM;
    scalar p_perp = sqrt(pE_CoM.x()*pE_CoM.x() + pE_CoM.y()*pE_CoM.y());
    if( p_perp > 1.e-10*mag_pE_CoM ) {
        scalar inv_p_perp = 1.0/p_perp;
        new_pE_CoM = vector( (pE_CoM.x() * pE_CoM.z() * sinCcosE - pE_CoM.y() * mag_pE_CoM * sinCsinE) * inv_p_perp + pE_CoM.x() * cosChi   ,\
                             (pE_CoM.y() * pE_CoM.z() * sinCcosE + pE_CoM.x() * mag_pE_CoM * sinCsinE) * inv_p_perp + pE_CoM.y() * cosChi   ,\
                              -p_perp * sinCcosE  +  pE_CoM.z() * cosChi);
    }
    else {
        new_pE_CoM = vector(mag_pE_CoM*sinCcosE,mag_pE_CoM*sinCsinE,mag_pE_CoM*cosChi);
    }

    vector new_pE = new_pE_CoM + vc_r*(gammac_SrqVcr*(vc_r&new_pE_CoM)+massE*gammaE_CoM*gamma_c);
    scalar new_gammaE = sqrt(1.0+sqr(mag(new_pE)/(massE*cu::c.value())));

    vector new_pN = -1.0 * new_pE_CoM + vc_r*(gammac_SrqVcr*(vc_r&(-1.0*new_pE_CoM))+massN*gammaN_CoM*gamma_c);
    scalar new_gammaN = sqrt(1.0+sqr(mag(new_pN)/(massN*cu::c.value())));

    Ue = new_pE/(massE*new_gammaE);
    Un = new_pN/(massN*new_gammaN);
}

// ************************************************************************* //
