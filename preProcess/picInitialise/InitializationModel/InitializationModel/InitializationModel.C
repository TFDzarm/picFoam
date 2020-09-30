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

#include "InitializationModel.H"
#include "fundamentalConstants.H"

template<class CloudType>
Foam::InitializationModel<CloudType>::InitializationModel(const dictionary& dict, CloudType& owner)
:
    dict_(dict),
    cloud_(owner),
    coeffDict_(dictionary::null),
    cellZoneIndex_(-1),
    numberDensities_(owner.typeIdList().size(),0.0),
    temperatures_(owner.typeIdList().size(),0.0)
{
}


template<class CloudType>
Foam::InitializationModel<CloudType>::InitializationModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    cloud_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    cellZoneIndex_(-1),
    numberDensities_(owner.typeIdList().size(),0.0),
    temperatures_(owner.typeIdList().size(),0.0)
{
    word cellZone = coeffDict_.lookupOrDefault<word>("cellZone","none");
    if(cellZone != "none")
        cellZoneIndex_ = owner.mesh().cellZones().findIndex(cellZone);
}


template<class CloudType>
Foam::InitializationModel<CloudType>::~InitializationModel()
{}

template<class CloudType>
void Foam::InitializationModel<CloudType>::readDefaultSettings()
{
    CloudType& cloud(this->cloud());


    if(!this->calculateNumberDensities())
    {
        const dictionary& numberDensitiesDict
        (
            dict_.subDict("numberDensities")
        );

        List<word> molecules(numberDensitiesDict.toc());

        forAll(molecules, i)
        {
            const word& moleculeName(molecules[i]);
            label typeId(findIndex(cloud.typeIdList(), moleculeName));
            if(typeId == -1)
                FatalErrorInFunction << "Species " << moleculeName << " is not defined in cloud" << nl << abort(FatalError);

            numberDensities_[typeId] = readScalar
            (
                numberDensitiesDict.lookup(molecules[i])
            );
        }
    }

    if(!this->calculateTemperatures())
    {
        const dictionary& temperatureDict
        (
            dict_.subDict("temperature")
        );

        List<word> moleculeTemperatures(temperatureDict.toc());

        bool electronVolt = dict_.lookupOrDefault<bool>("electronVolt",false);

        forAll(moleculeTemperatures, i)
        {
            const word& moleculeName(moleculeTemperatures[i]);
            label typeId(findIndex(cloud.typeIdList(), moleculeName));

            if(typeId == -1)
                FatalErrorInFunction << "Species " << moleculeName << " is not defined in cloud" << nl << abort(FatalError);

            temperatures_[typeId] = readScalar
            (
                temperatureDict.lookup(moleculeTemperatures[i])
            );

            if(electronVolt)
                temperatures_[typeId] *= constant::electromagnetic::e.value()/constant::physicoChemical::k.value();
        }
    }
}


template<class CloudType>
const CloudType&
Foam::InitializationModel<CloudType>::cloud() const
{
    return cloud_;
}


template<class CloudType>
CloudType&
Foam::InitializationModel<CloudType>::cloud()
{
    return cloud_;
}


template<class CloudType>
const Foam::dictionary&
Foam::InitializationModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::InitializationModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
Field<scalar>& Foam::InitializationModel<CloudType>::numberDensities()
{
    return numberDensities_;
}

template<class CloudType>
const Field<scalar>& Foam::InitializationModel<CloudType>::numberDensities() const
{
    return numberDensities_;
}

template<class CloudType>
Field<scalar>& Foam::InitializationModel<CloudType>::temperatures()
{
    return temperatures_;
}

template<class CloudType>
const Field<scalar>& Foam::InitializationModel<CloudType>::temperatures() const
{
    return temperatures_;
}

template<class CloudType>
const Foam::label& Foam::InitializationModel<CloudType>::cellZone() const
{
    return cellZoneIndex_;
}

#include "InitializationModelNew.C"
