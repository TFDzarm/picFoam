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

#include "DebyeLengthInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DebyeLengthInfo<CloudType>::DebyeLengthInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    vSqr_(cloud.typeIdList().size(),0.0),
    vDrift_(cloud.typeIdList().size(),Zero),
    nParticleTypes_(cloud.typeIdList().size(),0.0),
    sCharge_(cloud.typeIdList().size(),0.0),
    totalVolume_(sum(cloud.mesh().cellVolumes())),
    accountForDrift_(true),//default
    gCellLengthScale_(gAverage(cloud.cellLengthScale())),
    warnCellSize_(false)
{
    const polyMesh& mesh(cloud.mesh());
    totalVolume_ = sum(mesh.cellVolumes());
    reduce(totalVolume_,sumOp<scalar>());

    accountForDrift_.readIfPresent("accountForDrift",this->coeffDict());
    if(accountForDrift_)
        Info << "       |= Do not incorporate dirft velocities in temperature calculation" << endl;
    else
        Info << "       |= Incorporate the drift velocities in calculations in temperature calculation" << endl;

    warnCellSize_.readIfPresent("warnCellSize",this->coeffDict());
    if(warnCellSize_)
        Info << "       |= Warn when the average cell size is larger than the Debye length" << endl;
}

template<class CloudType>
Foam::DebyeLengthInfo<CloudType>::DebyeLengthInfo(const DebyeLengthInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    vSqr_(im.vSqr_),
    vDrift_(im.vDrift_),
    nParticleTypes_(im.nParticleTypes_),
    sCharge_(im.sCharge_),
    totalVolume_(im.totalVolume_),
    accountForDrift_(im.accountForDrift_),
    gCellLengthScale_(im.gCellLengthScale_),
    warnCellSize_(im.warnCellSize_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DebyeLengthInfo<CloudType>::~DebyeLengthInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DebyeLengthInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    //Calculate all required properties
    vSqr_[p.typeId()] += (p.U()& p.U())*p.nParticle();
    vDrift_[p.typeId()] += p.U()*p.nParticle();
    sCharge_[p.typeId()] += p.charge()*p.nParticle();
    nParticleTypes_[p.typeId()] += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::DebyeLengthInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());
    
    //Parallel COM the properties
    Pstream::listCombineGather(vSqr_, plusEqOp<scalar>());
    Pstream::listCombineGather(vDrift_, plusEqOp<vector>());
    Pstream::listCombineGather(sCharge_, plusEqOp<scalar>());
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());

    scalar nMoles = sum(nParticleTypes_);
    if(nMoles <= 0.0)
        return;

    scalar debye = 0.0;
    forAll(cloud.chargedSpecies(),i)
    {
        scalar id = cloud.chargedSpecies()[i];

        if(nParticleTypes_[id] <= 0.0)
            continue;

        scalar vDrift2 = 0;//default
        if(accountForDrift_)
            vDrift2 += (vDrift_[id]&vDrift_[id]);

        scalar T = cloud.constProps(id).mass()*(vSqr_[id] - vDrift2/nParticleTypes_[id]);
        T /= 3.0*constant::physicoChemical::k.value()*nParticleTypes_[id];

        scalar n = nParticleTypes_[id]/totalVolume_;
        scalar charge = sCharge_[id]/nParticleTypes_[id];
        debye += n*charge*charge/(constant::electromagnetic::epsilon0.value()*constant::physicoChemical::k.value()*T);
    }


    scalar debyeLength = 0.0;
    if(debye > 0.0)
        debyeLength = ::sqrt(1.0/debye);

    //Printing
    Info << "    Average Debye Length            = " << debyeLength << " m" << nl;

    if(warnCellSize_ && gCellLengthScale_ > debyeLength)
        Info << "    WARNING: Average cell size of " << gCellLengthScale_ << " m exceeds the Debye length!" << nl;

    //Reset list
    forAll(cloud.typeIdList(),i)
    {
        vSqr_[i] = 0.0;
        vDrift_[i] = Zero;
        nParticleTypes_[i] = 0.0;
        sCharge_[i] = 0.0;
    }
}

// ************************************************************************* //
