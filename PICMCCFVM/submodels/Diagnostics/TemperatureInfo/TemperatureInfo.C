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

#include "TemperatureInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TemperatureInfo<CloudType>::TemperatureInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    vSqr_(cloud.typeIdList().size(),0.0),
    vDrift_(cloud.typeIdList().size(),Zero),
    nParticleTypes_(cloud.typeIdList().size(),0.0),
    accountForDrift_(true),//default
    nMoles_(0.0)
{
    accountForDrift_.readIfPresent("accountForDrift",this->coeffDict());
}

template<class CloudType>
Foam::TemperatureInfo<CloudType>::TemperatureInfo(const TemperatureInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    vSqr_(im.vSqr_),
    vDrift_(im.vDrift_),
    nParticleTypes_(im.nParticleTypes_),
    accountForDrift_(im.accountForDrift_),
    nMoles_(im.nMoles_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TemperatureInfo<CloudType>::~TemperatureInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TemperatureInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    vSqr_[p.typeId()] += (p.U()& p.U())*p.nParticle();
    vDrift_[p.typeId()] += p.U()*p.nParticle();
    nParticleTypes_[p.typeId()] += p.nParticle();
    nMoles_ += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::TemperatureInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());
    reduce(nMoles_, sumOp<scalar>());

    scalar linearKineticEnergy = cloud.linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());
    Pstream::listCombineGather(vSqr_, plusEqOp<scalar>());
    Pstream::listCombineGather(vDrift_, plusEqOp<vector>());
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());

    Info << "    Average Temperature             = "
         << 2.0*linearKineticEnergy/(3.0*nMoles_*constant::physicoChemical::k.value()) << " K => " << 2.0*linearKineticEnergy/(3.0*nMoles_*constant::electromagnetic::e.value()) << " eV" << nl;
    forAll(cloud.typeIdList(),id)
    {
        if(nParticleTypes_[id] > 0)
        {
            vector vDrift = vDrift_[id]/nParticleTypes_[id];
            scalar vDrift2 = 0;
            if(accountForDrift_)
                vDrift2 = (vDrift&vDrift);
            scalar temperatur = cloud.constProps(id).mass()/(3.0*constant::physicoChemical::k.value())*(vSqr_[id]/nParticleTypes_[id] - vDrift2);
            Info << "    [" << cloud.typeIdList()[id] << "]                             = " << temperatur << " K == " << temperatur*constant::physicoChemical::k.value()/constant::electromagnetic::e.value() << " eV" << nl;

        }
    }

   //Reset list
   forAll(cloud.typeIdList(),i)
   {
       vSqr_[i] = 0.0;
       vDrift_[i] = Zero;
       nParticleTypes_[i] = 0.0;
   }
   nMoles_ = 0.0;

}

// ************************************************************************* //
