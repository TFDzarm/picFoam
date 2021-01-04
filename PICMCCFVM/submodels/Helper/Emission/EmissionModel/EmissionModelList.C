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
#include "EmissionModelList.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::EmissionModelList<CloudType>::EmissionModelList(CloudType& owner)
:
    PtrList<Foam::EmissionModel<CloudType>>()
{}


template<class CloudType>
Foam::EmissionModelList<CloudType>::EmissionModelList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<Foam::EmissionModel<CloudType>>()
{
    setupModels(dict, owner);//Setup the models
}

template<class CloudType>
void Foam::EmissionModelList<CloudType>::setupModels(const dictionary& dict, CloudType& owner)
{
    List<word> modelList(dict.lookup("EmissionModels"));

    //Check for multiple defintions so we construct only once
    forAll(modelList,i)
    {
        forAllReverse(modelList,j)
        {
            if(i==j)
                break;

            if(modelList[i] == modelList[j])
                FatalError << "EmissionModel " << modelList[i] << " is used mutiple times" << abort(FatalError);
        }

    }

    //Set the size of the model list
    this->setSize(modelList.size());

    //Construct all models
    forAll(modelList,i)
    {
        word modelType = modelList[i];

        this->set(
                  i,
                  EmissionModel<CloudType>::New(
                            modelType,
                            dict,
                            owner
                        )
                  );
    }
}

template<class CloudType>
Foam::EmissionModelList<CloudType>::EmissionModelList
(
    const EmissionModelList<CloudType>& iml
)
:
    PtrList<Foam::EmissionModel<CloudType>>(iml)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EmissionModelList<CloudType>::~EmissionModelList()
{}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::EmissionModelList<CloudType>::emission()
{
    //Call the emission function for ever model in the list
    forAll(*this, i)
        this->operator[](i).emission();
}
template<class CloudType>
void Foam::EmissionModelList<CloudType>::initilizeAll(label patchId, typename ParticleEmitter<CloudType>::VelocityModel model)
{
    //Call the initialisation function for every mode this will set up the emitter class
    forAll(*this, i)
        this->operator[](i).initilizeParticleEmitter(patchId,model);
}


template<class CloudType>
void Foam::EmissionModelList<CloudType>::collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    //Inform every model about the collision of particle p
    //Models will check if the particle is on an associated patch
    forAll(*this, i)
        this->operator[](i).collisionalEmission(p,td);
}
// ************************************************************************* //

