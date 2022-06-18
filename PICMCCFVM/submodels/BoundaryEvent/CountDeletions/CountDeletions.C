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

#include "CountDeletions.H"
#include "Pstream.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CountDeletions<CloudType>::CountDeletions
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryEvent<CloudType>(dict,cloud,associatedPatches),
    countParticle_(),
    countParcel_(),
    runningAverage_(),
    averageStart_(cloud.mesh().time().value())
{
    const polyMesh& mesh(cloud.mesh());

    countParticle_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
    countParcel_.setSize(mesh.boundaryMesh().size(),Field<label>(cloud.typeIdList().size(),0));
    runningAverage_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
}

template<class CloudType>
Foam::CountDeletions<CloudType>::CountDeletions(const CountDeletions<CloudType>& im)
:
    BoundaryEvent<CloudType>(im.owner_),
    countParticle_(im.countParticle_),
    runningAverage_(im.runningAverage_),
    averageStart_(im.averageStart_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CountDeletions<CloudType>::~CountDeletions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CountDeletions<CloudType>::collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    if(!td.keepParticle) {
        countParcel_[p.patch()][p.typeId()]++;
        countParticle_[p.patch()][p.typeId()]+= p.nParticle();
    }
}

template<class CloudType>
void Foam::CountDeletions<CloudType>::info()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    if(cloud.isInitializing())
        return;


    scalar dt = mesh.time().deltaTValue();
    scalar beta = dt/(mesh.time().value()-averageStart_);

    Info << nl << "[" << typeName << "]" << nl;

    forAll(this->associatedPatches(),i)
    {
        label patchId = this->associatedPatches()[i];
        const polyPatch& patch = mesh.boundaryMesh()[patchId];

        Info << "    " << patch.name() << ':' << nl;

        Pstream::listCombineGather(countParcel_[patchId], plusEqOp<label>());
        Pstream::listCombineScatter(countParcel_[patchId]);

        Pstream::listCombineGather(countParticle_[patchId], plusEqOp<scalar>());
        Pstream::listCombineScatter(countParticle_[patchId]);

        forAll(cloud.typeIdList(),j)
        {
            runningAverage_[patchId][j] = (1.0-beta)*runningAverage_[patchId][j] + beta*countParticle_[patchId][j];
            scalar rate = runningAverage_[patchId][j]/dt;

            Info << "        [" << cloud.typeIdList()[j] << "] parcel = " << countParcel_[patchId][j] << " | particle = " << countParticle_[patchId][j] << " | flow rate = " << rate << " 1/s" << nl;

            //Reset
            countParticle_[patchId][j] = 0.0;
            countParcel_[patchId][j] = 0;
        }
    }

}

// ************************************************************************* //
