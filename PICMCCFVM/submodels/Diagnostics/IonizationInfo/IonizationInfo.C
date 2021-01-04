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

#include "IonizationInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IonizationInfo<CloudType>::IonizationInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    ionizationStates_(cloud.typeIdList().size())
{
    //Check up to which limit we print info
    forAll(cloud.ionSpecies(), i)
    {
        label typeId = cloud.ionSpecies()[i];
        const dictionary& subDict = this->coeffDict().subDict(cloud.typeIdList()[typeId]);//Get a subdict with the species name

        label ionLimit = readLabel(subDict.lookup("ionizationLimit"));
        ionizationStates_[typeId].resize(ionLimit,0);
    }
}

template<class CloudType>
Foam::IonizationInfo<CloudType>::IonizationInfo(const IonizationInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    ionizationStates_(im.ionizationStates_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IonizationInfo<CloudType>::~IonizationInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IonizationInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    //If in the limit count the particle
    if(ionizationStates_[p.typeId()].size() >= p.Zstar())
    {
        ionizationStates_[p.typeId()][p.Zstar()-1]++;
    }
}

//- Print info
template<class CloudType>
void Foam::IonizationInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    forAll(cloud.ionSpecies(),i)
    {
        label typeId = cloud.ionSpecies()[i];
    
        //Parallel COM the list for the current species
        Pstream::listCombineGather(ionizationStates_[typeId],plusEqOp<label>());

        //Print the info
        word speciesName = cloud.typeIdList()[typeId];
        Info << "    Ionization states of " << speciesName << ":" << nl;

        forAll(ionizationStates_[typeId],stateI)
        {
            label count = ionizationStates_[typeId][stateI];
            if(count > 0)
                Info << "    " << stateI+1 << ": " << count << nl;
        }
        Info << nl;
    }

    //Reset list
    forAll(ionizationStates_,i) {
        if(ionizationStates_[i].size() > 0)
            ionizationStates_[i] = 0;
    }
}

// ************************************************************************* //

