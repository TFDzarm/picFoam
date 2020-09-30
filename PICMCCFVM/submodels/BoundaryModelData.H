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

Class
    Foam::BoundaryModelData

Description
    Helper class for BoundaryEvent and BoundaryModel classes

\*---------------------------------------------------------------------------*/

#ifndef BoundaryModelData_H
#define BoundaryModelData_H

#include "Istream.H"
#include "labelList.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class BoundaryModelData Declaration
\*---------------------------------------------------------------------------*/

// Forward declaration of classes
class BoundaryModelData;

// Forward declaration of friend functions
Istream& operator>>
(
    Istream& is,
    BoundaryModelData& pid
);


class BoundaryModelData
{
    // Private data

        //- Model name
        word modelName_;

        //- String including all patch names (patch1|patch2)
        word patchName_;

        //- List of all patchIds
        labelList patchIds_;

        //- Associated dictionary
        dictionary modelDict_;

public:

    // Constructor

        //- Construct null
        BoundaryModelData();


    // Member functions

        // Access

            //- Return const access to the reaction model name
            const word& modelName() const;

            //- Return const access to the patch name
            const word& patchName() const;

            //- Return const patchId list
            const labelList& patchIds() const;

            //- Return non-const patchId list
            labelList& patchIds();

            //- Return associated dictionary
            const dictionary& modelDict() const;


        // I-O

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                BoundaryModelData& pid
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

