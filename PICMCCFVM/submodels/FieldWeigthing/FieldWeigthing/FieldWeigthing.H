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
    Foam::FieldWeigthing

Description
    Base class for electric field weighting

SourceFiles
    FieldWeigthing.C

\*---------------------------------------------------------------------------*/

#ifndef FieldWeigthing_H
#define FieldWeigthing_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "barycentric.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tetIndices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FieldWeigthing Declaration
\*---------------------------------------------------------------------------*/

class FieldWeigthing
{
protected:

    //- Protected data

        //- The stored electric field
        volVectorField field_;

        //- Dictionary
        const dictionary& dict_;

        //- The mesh
        const fvMesh& mesh_;
public:

    //- Runtime type information
    TypeName("FieldWeigthing");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FieldWeigthing,
        dictionary,
        (
            const word& fieldName,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (fieldName, dict, mesh)
    );


    //- Constructors

        //- Construct from components
        FieldWeigthing
        (
            const word& fieldName,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Construct a copy
        FieldWeigthing(const FieldWeigthing& am);

        //- Construct and return a clone
        virtual autoPtr<FieldWeigthing> clone() const = 0;


    //- Selector
    static autoPtr<FieldWeigthing> New
    (
        const word& fieldName,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~FieldWeigthing();


    //- Return non-const the electric field
    volVectorField& field()
    {
        return field_;
    }

    //- Return const the electric field
    const volVectorField& field() const
    {
        return field_;
    }

    //- Return electric field vector for the coordinates
    virtual vector getFieldVector
    (
        const barycentric& coordinates,
        const tetIndices& tetIs
    ) const = 0;

    //- Update the interpolation
    virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#define makeFieldWeigthing()                                    \
                                                                               \
    defineTypeNameAndDebug                                        \
    (                                                                          \
        FieldWeigthing,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineRunTimeSelectionTable                                        \
    (                                                                          \
        FieldWeigthing,                                       \
        dictionary                                                             \
    );


#define makeFieldWeigthingType(SS)                            \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                     \
                                                                               \
    FieldWeigthing::                                          \
        adddictionaryConstructorToTable<SS>                         \
            add##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

