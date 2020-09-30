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

#include "KineticEnergyInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KineticEnergyInfo<CloudType>::KineticEnergyInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    kineticEnergyofTypes_(cloud.typeIdList().size(),0.0),
    relativistic_(false)
{
    relativistic_.readIfPresent("relativistic",this->coeffDict());

    if(relativistic_)
        Info << "       |= Calculated relativistically" << endl;
    else
        Info << "       |= Calculated classical" << endl;
}

template<class CloudType>
Foam::KineticEnergyInfo<CloudType>::KineticEnergyInfo(const KineticEnergyInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    kineticEnergyofTypes_(im.kineticEnergyofTypes_),
    relativistic_(im.relativistic_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KineticEnergyInfo<CloudType>::~KineticEnergyInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KineticEnergyInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const CloudType& cloud(this->owner());
    const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(p.typeId());

    if(relativistic_) {
        scalar gamma = 1.0/sqrt(1.0-sqr(mag(p.U())/constant::universal::c.value()));
        kineticEnergyofTypes_[p.typeId()] += (gamma-1.0)*cP.mass()*constant::universal::c.value()*constant::universal::c.value();
    }
    else {
        kineticEnergyofTypes_[p.typeId()] += 0.5*cP.mass()*(p.U() & p.U())*p.nParticle();
    }
}

//- Print info
template<class CloudType>
void Foam::KineticEnergyInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    scalar linearKineticEnergy = cloud.linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());
    Pstream::listCombineGather(kineticEnergyofTypes_, plusEqOp<scalar>());

    Info << "    Total linear kinetic energy     = "
         << linearKineticEnergy << " J == " << linearKineticEnergy/constant::electromagnetic::e.value() << " eV" << nl;
    forAll(cloud.typeIdList(),id)
    {
        scalar kineticEnergy = kineticEnergyofTypes_[id];
        Info << "    [" << cloud.typeIdList()[id] << "]                             = " << kineticEnergy << " J == " << kineticEnergy/constant::electromagnetic::e.value() << " eV" << nl;
    }
   //Reset list
   forAll(cloud.typeIdList(),i)
   {
       kineticEnergyofTypes_[i] = 0.0;
   }

}

// ************************************************************************* //
