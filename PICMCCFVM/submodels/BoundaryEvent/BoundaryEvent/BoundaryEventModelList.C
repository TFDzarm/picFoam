/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
#include <algorithm>
#include "BoundaryEventModelList.H"
#include "polyMesh.H"
#include "stringListOps.H"
#include "emptyPolyPatch.H"
#include "BoundaryModelData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::BoundaryEventModelList<CloudType>::BoundaryEventModelList(CloudType& owner)
:
    PtrList<Foam::BoundaryEvent<CloudType>>()
{}


template<class CloudType>
Foam::BoundaryEventModelList<CloudType>::BoundaryEventModelList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<Foam::BoundaryEvent<CloudType>>()
{
    setupModels(dict, owner);
}

/*
Foam::BoundaryEventModelList<CloudType>::setupModels

Read the modelList from constant/picProperties and construct the model
*/
template<class CloudType>
void Foam::BoundaryEventModelList<CloudType>::setupModels(const dictionary& dict, CloudType& owner)
{
    const polyMesh& mesh(owner.mesh());
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const wordList allPatchNames = bMesh.names();

    List<BoundaryModelData> modelData(dict.lookup("modelList"));
    if(modelData.empty())
        return;

    DynamicList<word> selectedModels;
    forAllReverse(modelData, i)
    {
        const word& patchName = modelData[i].patchName();
        labelList patchIDs = findStrings(patchName, allPatchNames);
        label patchCount = std::count(patchName.begin(),patchName.end(),'|') + 1;

        if (patchIDs.empty() || patchIDs.size() < patchCount )
        {
            WarningInFunction
                << "Cannot find one or more patch names matching (" << patchName << ")"
                << endl;
        }
        modelData[i].patchIds().transfer(patchIDs);

        label index = findIndex(selectedModels,modelData[i].modelName());
        if(index == -1)
            selectedModels.append(modelData[i].modelName());
    }

    this->setSize(selectedModels.size());//Set the size of the pointer list
    forAll(selectedModels,i)
    {
        word modelType = selectedModels[i];

        DynamicList<label> associatedPatches;//List of patches this model is defined on, will be saved by the model
        forAll(modelData,j)
        {
            if(modelData[j].modelName() == modelType)
                associatedPatches.append(modelData[j].patchIds());
        }

        //Construct the model
        this->set(i,
                  BoundaryEvent<CloudType>::New(
                      modelType,
                      dict,
                      owner,
                      associatedPatches
                      ));
    }
}

template<class CloudType>
Foam::BoundaryEventModelList<CloudType>::BoundaryEventModelList
(
    const BoundaryEventModelList<CloudType>& iml
)
:
    PtrList<Foam::BoundaryEvent<CloudType>>(iml)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundaryEventModelList<CloudType>::~BoundaryEventModelList()
{}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
Foam::BoundaryEventModelList<CloudType>::onCollision

A particle has collided with a patch, check if we have a model on this patch
*/
template<class CloudType>
void Foam::BoundaryEventModelList<CloudType>::onCollision
(
    typename CloudType::parcelType& p,
    typename CloudType::parcelType::trackingData& td
)
{
    forAll(*this, i)
    {
        if(this->operator[](i).interactWithPatch(p.patch()))//Should we interact with this patch, checks associatedPatches
           this->operator[](i).collisionEvent(p, td);//Interact...
    }
}

/*
Foam::BoundaryEventModelList<CloudType>::onCollision

A particle was ejected from a patch, check if we have a model on this patch
*/
template<class CloudType>
void Foam::BoundaryEventModelList<CloudType>::onEjection
(
    typename CloudType::parcelType& p,
    label patchId
)
{
    forAll(*this, i)
    {
        if(this->operator[](i).interactWithPatch(patchId))//Should we interact with this patch, checks associatedPatches
            this->operator[](i).ejectionEvent(p, patchId);//Interact...
    }
}

/*
Foam::BoundaryEventModelList<CloudType>::info

The timestep has finished, print some info...
*/
template<class CloudType>
void Foam::BoundaryEventModelList<CloudType>::info()
{
    forAll(*this, i)
    {
        this->operator[](i).info();
    }
}

/*
Foam::BoundaryEventModelList<CloudType>::postMove

All parcels have been moved, do something...
*/
template<class CloudType>
void Foam::BoundaryEventModelList<CloudType>::postMove()
{
    forAll(*this, i)
    {
        this->operator[](i).postMove();
    }
}

// ************************************************************************* //

