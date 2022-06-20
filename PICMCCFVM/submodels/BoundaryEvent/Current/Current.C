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

#include "Current.H"
#include "Pstream.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Current<CloudType>::Current
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryEvent<CloudType>(dict,cloud,associatedPatches),
    chargeCounter_(),
    runningAverage_(),
    averageStart_(cloud.mesh().time().value())
{
    const polyMesh& mesh(cloud.mesh());

    chargeCounter_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
    runningAverage_.setSize(mesh.boundaryMesh().size(),Field<scalar>(cloud.typeIdList().size(),0.0));
}

template<class CloudType>
Foam::Current<CloudType>::Current(const Current<CloudType>& im)
:
    BoundaryEvent<CloudType>(im.owner_),
    chargeCounter_(im.chargeCounter_),
    runningAverage_(im.runningAverage_),
    averageStart_(im.averageStart_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Current<CloudType>::~Current()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::Current<CloudType>::collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    //Charge leaving the system
    if(!td.keepParticle) {
        chargeCounter_[p.patch()][p.typeId()]+= p.nParticle()*p.charge();
    }
}

template<class CloudType>
void Foam::Current<CloudType>::ejectionEvent(typename CloudType::parcelType& p, label patchId)
{
    //Charge entering the system
    chargeCounter_[patchId][p.typeId()]-= p.nParticle()*p.charge();
}

template<class CloudType>
void Foam::Current<CloudType>::info()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    //Do not operate when initializing, time will be 0
    if(cloud.isInitializing())
        return;


    scalar dt = mesh.time().deltaTValue();
    scalar beta = dt/(mesh.time().value()-averageStart_);

    Info << nl << "[" << typeName << "]\n";
    forAll(this->associatedPatches(),i)
    {
        //Get patch
        label patchId = this->associatedPatches()[i];
        const polyPatch& patch = mesh.boundaryMesh()[patchId];

        Info << "    " << patch.name() << ':' << nl;

        //Parallel communication
        Pstream::listCombineGather(chargeCounter_[patchId], plusEqOp<scalar>());
        Pstream::listCombineScatter(chargeCounter_[patchId]);

        scalar totalCurrent = 0.0;
        forAll(cloud.chargedSpecies(),j)//Only consider charged species
        {
            //Get the actual typeId from chargedSpecies list
            label typeId = cloud.chargedSpecies()[j];

            //Calculate running average of the charge moving across the patch(patchId)
            runningAverage_[patchId][typeId] = (1.0-beta)*runningAverage_[patchId][typeId] + beta*chargeCounter_[patchId][typeId];

            scalar rate = runningAverage_[patchId][typeId]/dt;
            //Total current over all species
            totalCurrent +=  rate;
            Info << "        [" << cloud.typeIdList()[typeId] << "]: " << rate << " A" << nl;

            //Reset
            chargeCounter_[patchId][typeId] = 0.0;
        }
        Info << "        Total current = " << totalCurrent << " A" << nl;
    }

}

// ************************************************************************* //
