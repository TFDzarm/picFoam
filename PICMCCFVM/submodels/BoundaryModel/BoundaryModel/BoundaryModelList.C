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
#include <algorithm>
#include "BoundaryModelList.H"
#include "polyMesh.H"
#include "stringListOps.H"
#include "emptyPolyPatch.H"
#include "BoundaryModelData.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::BoundaryModelList<CloudType>::BoundaryModelList(CloudType& owner)
:
    PtrList<Foam::BoundaryModel<CloudType>>()
{}


template<class CloudType>
Foam::BoundaryModelList<CloudType>::BoundaryModelList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<Foam::BoundaryModel<CloudType>>()
{
    setupModels(dict,owner);
}

template<class CloudType>
void Foam::BoundaryModelList<CloudType>::setupModels(const dictionary& dict, CloudType& owner)
{
    const polyMesh& mesh(owner.mesh());
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const wordList allPatchNames = bMesh.names();

    //Get number of patches
    DynamicList<label> polyPatches;
    forAll(owner.mesh().boundaryMesh(), patchId)
    {
        const polyPatch& patch = owner.mesh().boundaryMesh()[patchId];
        if(isType<polyPatch>(patch))
            polyPatches.append(patchId);
    }


    List<BoundaryModelData> modelData(dict.lookup("modelList"));
    //Check if boundarymodels are definied
    if(modelData.empty() && polyPatches.size() > 0)
    {
        FatalError << "No BoundaryModel was provided for patches" << nl;
                      forAll(polyPatches,i)
                              FatalError << "    " << allPatchNames[polyPatches[i]] << nl;
        FatalError << abort(FatalError);
    }

    //Get patchIds and check if no boundarymodel is used multiple times
    DynamicList<word> selectedModels;
    forAllReverse(modelData, i)
    {
        const word& patchName = modelData[i].patchName();

        labelList patchIDs = findStrings(patchName, allPatchNames);
        label patchCount = std::count(patchName.begin(),patchName.end(),'|') + 1;

        if (patchIDs.empty() || patchIDs.size() < patchCount )
        {
            FatalError
                << "Cannot find one or more patch names matching (" << patchName << ")"
                << abort(FatalError);
        }
        modelData[i].patchIds().transfer(patchIDs);


        label index = findIndex(selectedModels,modelData[i].modelName());

        if(index == -1)
            selectedModels.append(modelData[i].modelName());
        else
        {
            //Allow same model on different patches only do the check for modelType none
            if(modelData[i].modelName() == "none") {
                FatalError << "Found separate BoundaryModel definitions for model " << modelData[i].modelName() << nl
                           << abort(FatalError);
            }
        }
    }

    //Check if polypatch was definied multiple times
    forAll(modelData,i)
    {
        forAll(modelData,j)
        {
            if(i==j)
                break;
            forAll(modelData[j].patchIds(),idJ)
            {
                label index = findIndex(modelData[i].patchIds(),modelData[j].patchIds()[idJ]);
                if(index != -1)
                {
                    FatalError << "Patch " << allPatchNames[modelData[i].patchIds()[index]] << " is used in mutiple BoundaryModels" << abort(FatalError);
                }
            }

        }

    }

    //Check if all polyPatches have a boundarymodel
    forAll(polyPatches,i)
    {
        label patchId = polyPatches[i];
        bool foundEntry = false;
        forAll(modelData,i)
        {
            foundEntry = findIndex(modelData[i].patchIds(),patchId) != -1 ? true : false;
            if(foundEntry)
                break;
        }
        if(!foundEntry)
        {
            FatalError << "No BoundaryModel definition for patch " << allPatchNames[patchId] << abort(FatalError);
        }
    }

    this->setSize(modelData.size());
    forAll(modelData,i)
    {
        word modelType = modelData[i].modelName();

        DynamicList<label> associatedPatches;

        associatedPatches.append(modelData[i].patchIds());

        this->set(i,
                  BoundaryModel<CloudType>::New(
                      modelType,
                      modelData[i].modelDict(),
                      owner,
                      associatedPatches
                      ));
    }
}

template<class CloudType>
Foam::BoundaryModelList<CloudType>::BoundaryModelList
(
    const BoundaryModelList<CloudType>& iml
)
:
    PtrList<Foam::BoundaryModel<CloudType>>(iml)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundaryModelList<CloudType>::~BoundaryModelList()
{}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BoundaryModelList<CloudType>::injection()
{
    forAll(*this, i)
    {
        this->operator[](i).injection();
    }
}

template<class CloudType>
void Foam::BoundaryModelList<CloudType>::preUpdate_Boundary()
{
    forAll(*this, i)
    {
        this->operator[](i).preUpdate_Boundary();
    }
}

template<class CloudType>
void Foam::BoundaryModelList<CloudType>::postUpdate_Boundary()
{
    forAll(*this, i)
    {
        this->operator[](i).postUpdate_Boundary();
    }
}

template<class CloudType>
bool Foam::BoundaryModelList<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    forAll(*this, i)
    {
        if(this->operator[](i).interactWithPatch(p.patch()))
        {
            return this->operator[](i).particleBC(p,td);//Only one patch will be valid
        }
    }
    return false;//default
}

template<class CloudType>
void Foam::BoundaryModelList<CloudType>::onEjection(typename CloudType::parcelType& p, label patchId)
{
    forAll(*this, i)
    {
        if(this->operator[](i).interactWithPatch(patchId))
        {
            this->operator[](i).particleEjection(p, patchId);
        }
    }
}



// ************************************************************************* //

