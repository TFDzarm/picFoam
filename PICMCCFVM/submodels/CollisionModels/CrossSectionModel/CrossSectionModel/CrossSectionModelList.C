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
#include "CrossSectionModelList.H"
#include "polyMesh.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionList<CloudType,Type>::CrossSectionList(CloudType& owner)
:
    PtrList<Foam::CrossSectionModel<CloudType,Type>>()
{}

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionList<CloudType,Type>::CrossSectionList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<Foam::CrossSectionModel<CloudType,Type>>()
{
    //There can be one model for every species
    this->setSize(owner.typeIdList().size());

    //Lookup which models were specified by the user
    forAllConstIter(IDLList<entry>, dict, iter)
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();
            label typeId = findIndex(owner.typeIdList(),iter().keyword());
            if(typeId == -1)
            {
                FatalErrorInFunction << "No typeId named " << iter().keyword() << abort(FatalError);
            }
            if(this->set(typeId))
            {
                FatalErrorInFunction << "Model for " << iter().keyword() << " is already defined" << abort(FatalError);
            }
            //Based on the template "Type" read the correct entry
            word modelName;
            if(Type == Foam::crossSectionType::ElectronElasticCS)
                modelName = subDict.lookup<word>("ElasticCrossSection");
            else if(Type == Foam::crossSectionType::ElectronExciationCS)
                modelName = subDict.lookup<word>("ExcitationCrossSection");
            else if(Type == Foam::crossSectionType::ElectronIonizationCS)
                modelName = subDict.lookup<word>("IonizationCrossSection");
            //Construct the model
            this->set(
                          typeId,
                          CrossSectionModel<CloudType,Type>::New(
                              modelName,
                              subDict,
                              owner,
                              typeId
                          )
                      );
        }

    }
}



template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionList<CloudType,Type>::CrossSectionList
(
    const CrossSectionList<CloudType,Type>& iml
)
:
    PtrList<Foam::CrossSectionModel<CloudType,Type>>(iml)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::CrossSectionList<CloudType,Type>::~CrossSectionList()
{}




// ************************************************************************* //

