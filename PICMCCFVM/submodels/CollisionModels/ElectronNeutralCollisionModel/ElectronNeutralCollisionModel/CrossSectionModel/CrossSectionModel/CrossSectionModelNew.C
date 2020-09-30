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

#include "CrossSectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType, Foam::crossSectionType Type>
Foam::autoPtr<Foam::CrossSectionModel<CloudType,Type>>
Foam::CrossSectionModel<CloudType,Type>::New
(
    const word& modelName,
    const dictionary& dict,
    CloudType& owner,
    const label associatedTypeId
)
{
    word type;
    if(Type ==  Foam::crossSectionType::ElectronElasticCS)
        type = "Elastic";
    else if(Type ==  Foam::crossSectionType::ElectronExciationCS)
        type = "Exciation";
    else if(Type ==  Foam::crossSectionType::ElectronIonizationCS)
        type = "Ionization";

    Info<< "|->    Selecting " << type << "-CrossSectionModel " << modelName << " for typeid " << dict.dictName() << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown CrossSectionModel type "
            << modelName << nl << nl
            << "Valid CrossSectionModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<CrossSectionModel<CloudType,Type>>(cstrIter()(dict, owner, associatedTypeId));
}

// ************************************************************************* //
