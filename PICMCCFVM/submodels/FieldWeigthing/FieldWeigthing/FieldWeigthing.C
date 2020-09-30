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

#include "FieldWeigthing.H"
#include "runTimeSelectionTables.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FieldWeigthing::FieldWeigthing
(
    const word& fieldName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    field_(
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh    
    ),
    dict_(dict),
    mesh_(mesh)
{}


Foam::FieldWeigthing::FieldWeigthing
(
    const FieldWeigthing& am
)
:
    field_(am.field_),
    dict_(am.dict_),
    mesh_(am.mesh_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::FieldWeigthing>
Foam::FieldWeigthing::New
(
    const word& fieldName,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word type(dict.lookup("FieldWeigthing"));

    Info<< "+ Selecting FieldWeigthing model " << type << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown FieldWeigthing model " << type
            << ", constructor not in hash table" << nl << nl
            << "    Valid FieldWeigthing models are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<FieldWeigthing>(cstrIter()(fieldName, dict, mesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FieldWeigthing::~FieldWeigthing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FieldWeigthing::update()
{}

// ************************************************************************* //
