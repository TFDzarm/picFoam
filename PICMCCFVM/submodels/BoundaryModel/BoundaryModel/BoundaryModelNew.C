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

#include "BoundaryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
Foam::BoundaryModel<CloudType>::New Called in setupModels of BoundaryModelList
Create and return the pointer
*/
template<class CloudType>
Foam::autoPtr<Foam::BoundaryModel<CloudType>>
Foam::BoundaryModel<CloudType>::New
(
        const word& modelType,
        const dictionary& dict,
        CloudType& owner,
        const List<label>& associatedPatches
)
{
    const polyMesh& mesh(owner.mesh());
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const wordList allPatchNames = bMesh.names();

    Info<< "+ Selecting BoundaryModel " << modelType << " on patches: " << nl;
    forAll(associatedPatches,i)
            Info << "|->    " << allPatchNames[associatedPatches[i]] << nl;
    Info << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    //Model does not exist, print all known models
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown BoundaryModel type "
            << modelType << nl << nl
            << "Valid BoundaryModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<BoundaryModel<CloudType>>(cstrIter()(dict, owner, associatedPatches));
}


// ************************************************************************* //
