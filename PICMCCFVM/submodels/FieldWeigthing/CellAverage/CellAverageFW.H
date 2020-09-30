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
    Foam::CellAverageFW

Description
    Basic averaging model. Electric field is constant across a cell.

SourceFiles
    CellAverageFW.C

\*---------------------------------------------------------------------------*/

#ifndef CellAverageFW_H
#define CellAverageFW_H

#include "FieldWeigthing.H"
#include "pointMesh.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class CellAverageFW Declaration
\*---------------------------------------------------------------------------*/

class CellAverageFW
:
    public FieldWeigthing
{
public:

    //- Runtime type information
    TypeName("CellAverage");


    //- Constructors

        //- Construct from components
        CellAverageFW
        (
            const word& fieldName,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        CellAverageFW(const CellAverageFW& am);

        //- Construct and return a clone
        virtual autoPtr<FieldWeigthing> clone() const
        {
            return autoPtr<FieldWeigthing>
            (
                new CellAverageFW(*this)
            );
        }


    //- Destructor
    virtual ~CellAverageFW();


    //- Member Functions

        //- Return electric field vector for the coordinates
        vector getFieldVector
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
