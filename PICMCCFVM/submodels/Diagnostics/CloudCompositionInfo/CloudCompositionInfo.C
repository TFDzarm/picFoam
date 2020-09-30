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

#include "CloudCompositionInfo.H"
#include "constants.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudCompositionInfo<CloudType>::CloudCompositionInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    nParcelsTypes_(cloud.typeIdList().size(),0),
    nParticleTypes_(cloud.typeIdList().size(),0.0),
    totalVolume_(0.0),
    printNumberDensity_(false)
{
    printNumberDensity_.readIfPresent("printNumberDensity",this->coeffDict());

    const polyMesh& mesh(cloud.mesh());
    totalVolume_ = sum(mesh.cellVolumes());
    reduce(totalVolume_,sumOp<scalar>());
}

template<class CloudType>
Foam::CloudCompositionInfo<CloudType>::CloudCompositionInfo(const CloudCompositionInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    nParcelsTypes_(im.nParcelsTypes_),
    nParticleTypes_(im.nParticleTypes_),
    totalVolume_(im.totalVolume_),
    printNumberDensity_(im.printNumberDensity_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudCompositionInfo<CloudType>::~CloudCompositionInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudCompositionInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    nParcelsTypes_[p.typeId()]++;
    nParticleTypes_[p.typeId()] += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::CloudCompositionInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());
    Pstream::listCombineGather(nParcelsTypes_,plusEqOp<label>());
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());

    scalar nParticle = sum(nParticleTypes_);
    Info << "    Number of molecules             = "
         << nParticle
         << nl;
    forAll(cloud.typeIdList(),id)
    {
        if(nParcelsTypes_[id]) {
            Info<< "    [" << cloud.typeIdList()[id] << "]" << " parcels | particles        = " << nParcelsTypes_[id] << " | " << nParticleTypes_[id];
            if(printNumberDensity_)
                Info << " -> " << nParticleTypes_[id]/totalVolume_ << " m^-3" << nl;
            else
                Info << nl;
        }
    }


    //Reset list
    forAll(nParcelsTypes_,i) {
        nParcelsTypes_[i] = 0;
        nParticleTypes_[i] = 0.0;
    }
}

// ************************************************************************* //
