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
    totalCountParticle_(),
    countParticle_(),
    runningAverage_(),
    averageStart_(cloud.mesh().time().value())
{
    const polyMesh& mesh(cloud.mesh());

    countParticle_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
    totalCountParticle_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
    runningAverage_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
}

template<class CloudType>
Foam::CountDeletions<CloudType>::CountDeletions(const CountDeletions<CloudType>& im)
:
    BoundaryEvent<CloudType>(im.owner_),
    totalCountParticle_(im.totalCountParticle_),
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
    //Count the number of deleted particles
    if(!td.keepParticle) {
        countParticle_[p.patch()][p.typeId()]+= p.nParticle();
    }
}

template<class CloudType>
void Foam::CountDeletions<CloudType>::info()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    if(cloud.isInitializing())//Do nothing when we are initializing
        return;


    //Calculate a running average...
    scalar dt = mesh.time().deltaTValue();
    scalar beta = dt/(mesh.time().value()-averageStart_);

    forAll(this->associatedPatches(),i)
    {
        label patchId = this->associatedPatches()[i];
        const polyPatch& patch = mesh.boundaryMesh()[patchId];

        //Parallel COM
        Pstream::listCombineGather(countParticle_[patchId], plusEqOp<scalar>());
        Pstream::listCombineScatter(countParticle_[patchId]);

        forAll(cloud.typeIdList(),j)
        {
            const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(j);
            totalCountParticle_[patchId][j] += countParticle_[patchId][j];

            runningAverage_[patchId][j] = (1.0-beta)*runningAverage_[patchId][j] + beta*countParticle_[patchId][j];//running average

            scalar rate = runningAverage_[patchId][j]/dt;
            Info << "[patch: " << patch.name() << "]" << " Deletion rate (" << cloud.typeIdList()[j] << "): " << rate << " particle/s; I: " << rate * cP.charge() << " A" << endl;

            scalar total = totalCountParticle_[patchId][j];
            Info << "[patch: " << patch.name() << "]" << " Total deletion (" << cloud.typeIdList()[j] << "): " << total << " particles" << endl;

            countParticle_[patchId][j] = 0.0;
        }
    }

}

// ************************************************************************* //
