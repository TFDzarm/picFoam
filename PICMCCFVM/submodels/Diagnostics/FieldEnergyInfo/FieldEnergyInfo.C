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

#include "FieldEnergyInfo.H"
#include "volFieldsFwd.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FieldEnergyInfo<CloudType>::FieldEnergyInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    fieldEnergyofTypes_(cloud.typeIdList().size(),0.0),
    fromParcelDistribution_(false)
{
    fromParcelDistribution_.readIfPresent("fromDistribution",this->coeffDict());

    if(fromParcelDistribution_)
        Info << "       |= Calculated from parcel distribution" << endl;
    else
        Info << "       |= Calculated from electric field" << endl;
}

template<class CloudType>
Foam::FieldEnergyInfo<CloudType>::FieldEnergyInfo(const FieldEnergyInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    fieldEnergyofTypes_(im.fieldEnergyofTypes_),
    fromParcelDistribution_(im.fromParcelDistribution_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FieldEnergyInfo<CloudType>::~FieldEnergyInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FieldEnergyInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const CloudType& cloud(this->owner());
    const volScalarField& phiE(cloud.elpotentialField());

    if(fromParcelDistribution_)
        fieldEnergyofTypes_[p.typeId()] += p.charge()*p.nParticle()*phiE[p.cell()];
}

//- Print info
template<class CloudType>
void Foam::FieldEnergyInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    scalar fieldEnergy = 0;

    if(fromParcelDistribution_)
    {
        const volScalarField& phiE(cloud.elpotentialField());
        const volScalarField& spaceChargeDensity(cloud.spaceChargeDensity());

        scalar spaceChargeEnergy = 0.0;

        forAll(cloud.mesh().cells(),celli)
            spaceChargeEnergy+= phiE[celli]*spaceChargeDensity[celli]*cloud.mesh().cellVolumes()[celli];

        reduce(spaceChargeEnergy, sumOp<scalar>());
        Pstream::listCombineGather(fieldEnergyofTypes_, plusEqOp<scalar>());
        fieldEnergy = 0.5*sum(fieldEnergyofTypes_) + spaceChargeEnergy;

    }
    else
    {
        const volVectorField& Efield(cloud.electricField());

        forAll(cloud.mesh().cells(),celli)
            fieldEnergy += magSqr(Efield[celli])*cloud.mesh().cellVolumes()[celli];

        reduce(fieldEnergy, sumOp<scalar>());
        fieldEnergy *= 0.5*constant::electromagnetic::epsilon0.value();
    }

    Info << "    Total field energy              = "
         << fieldEnergy << " J == " << fieldEnergy/constant::electromagnetic::e.value() << " eV" << nl;
    if(fromParcelDistribution_)
    {
        forAll(cloud.typeIdList(),id)
        {
            scalar typeFieldEnergy = fieldEnergyofTypes_[id];
            Info << "    [" << cloud.typeIdList()[id] << "]                             = " << typeFieldEnergy << " J == " << typeFieldEnergy/constant::electromagnetic::e.value() << " eV" << nl;
        }

        //Reset list
        forAll(cloud.typeIdList(),i)
            fieldEnergyofTypes_[i] = 0.0;
    }
}

// ************************************************************************* //
