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

#include "NanbuElectronNeutralCollision.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuElectronNeutralCollision<CloudType>::NanbuElectronNeutralCollision
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ElectronNeutralCollisionModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanbuElectronNeutralCollision<CloudType>::~NanbuElectronNeutralCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NanbuElectronNeutralCollision<CloudType>::active() const
{
    return true;
}

template<class CloudType>
void Foam::NanbuElectronNeutralCollision<CloudType>::updateVelocity(scalar eV, typename CloudType::parcelType& pE, typename CloudType::parcelType& pN)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    scalar massE = cloud.constProps(pE.typeId()).mass();
    scalar massN = cloud.constProps(pN.typeId()).mass();

    scalar phi = constant::mathematical::twoPi*rndGen.scalar01();
    scalar cosChi;
    if(eV < 1e-30)
        cosChi = 1.0;
    else
        cosChi = 1.0/eV*(eV+2.-2.*::pow(1.+eV,rndGen.scalar01()));//Scattering angle

    scalar sinChi;
    if(cosChi*cosChi > 1.0)
    {
        cosChi = 1.0;
        sinChi = 0.0;
    }
    else
        sinChi = sqrt( 1.0 - cosChi*cosChi );

    //Relative velocity
    vector g = pE.U() - pN.U();
    scalar g_mag = mag(g);
    scalar g_perp = ::sqrt(::pow(g.y(),2)+::pow(g.z(),2));

    //Update the velocity according to Nanbu
    vector h = vector(g_perp*::cos(phi) , -1.0*(g.y()*g.x()*::cos(phi)+g_mag*g.z()*::sin(phi))/g_perp , -1.0*(g.z()*g.x()*::cos(phi)-g_mag*g.y()*::sin(phi))/g_perp);

    pE.U() = pE.U() - massN/(massN+massE)*(g*(1.0-cosChi)+h*sinChi);
    pN.U() = pN.U() + massE/(massN+massE)*(g*(1.0-cosChi)+h*sinChi);
}

//Same as above, but uses sampled velocity Un to calculate the scattering
template<class CloudType>
void Foam::NanbuElectronNeutralCollision<CloudType>::updateVelocity(scalar eV, typename CloudType::parcelType& pE, vector& Un, label idN)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    scalar massE = cloud.constProps(pE.typeId()).mass();
    scalar massN = cloud.constProps(idN).mass();

    scalar phi = constant::mathematical::twoPi*rndGen.scalar01();
    scalar cosChi;
    if(eV < 1e-30)
        cosChi = 1.0;
    else
        cosChi = 1.0/eV*(eV+2.-2.*::pow(1.+eV,rndGen.scalar01()));

    scalar sinChi;
    if(cosChi*cosChi > 1.0)
    {
        cosChi = 1.0;
        sinChi = 0.0;
    }
    else
        sinChi = sqrt( 1.0 - cosChi*cosChi );

    vector g = pE.U() - Un;
    scalar g_mag = mag(g);
    scalar g_perp = ::sqrt(::pow(g.y(),2)+::pow(g.z(),2));

    vector h = vector(g_perp*::cos(phi) , -1.0*(g.y()*g.x()*::cos(phi)+g_mag*g.z()*::sin(phi))/g_perp , -1.0*(g.z()*g.x()*::cos(phi)-g_mag*g.y()*::sin(phi))/g_perp);

    pE.U() = pE.U() - massN/(massN+massE)*(g*(1.0-cosChi)+h*sinChi);
    Un = Un + massE/(massN+massE)*(g*(1.0-cosChi)+h*sinChi);
}

// ************************************************************************* //
