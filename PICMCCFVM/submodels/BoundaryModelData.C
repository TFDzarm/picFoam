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

#include "BoundaryModelData.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::BoundaryModelData::BoundaryModelData()
:
    modelName_("unknownModelName"),
    patchName_("unknownPatch"),
    patchIds_(),
    modelDict_(dictionary::null)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::BoundaryModelData::modelName() const
{
    return modelName_;
}


const Foam::word& Foam::BoundaryModelData::patchName() const
{
    return patchName_;
}

const Foam::labelList& Foam::BoundaryModelData::patchIds() const
{
    return patchIds_;
}

Foam::labelList& Foam::BoundaryModelData::patchIds()
{
    return patchIds_;
}

const Foam::dictionary& Foam::BoundaryModelData::modelDict() const
{
    return modelDict_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    BoundaryModelData& pid
)
{
    is.check("Istream& operator>>(Istream&, patchInteractionData&)");

    const dictionaryEntry entry(dictionary::null, is);

    pid.modelDict_ = entry.dict();
    pid.patchName_ = entry.keyword();
    entry.lookup("type") >> pid.modelName_;

    return is;
}


// ************************************************************************* //

