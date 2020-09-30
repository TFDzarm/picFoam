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
    Foam::VolumeWeightingFW

Description
    Use Foam::interpolationCellPoint to do a linear interpolation between cell - vertex data and the particle position

SourceFiles
    VolumeWeightingFW.C

\*---------------------------------------------------------------------------*/

#ifndef VolumeWeightingFW_H
#define VolumeWeightingFW_H

#include "FieldWeigthing.H"
#include "autoPtr.H"
#include "pointMesh.H"
#include "tetIndices.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class Basic Declaration
\*---------------------------------------------------------------------------*/

class VolumeWeightingFW
:
    public FieldWeigthing
{
    //- Pointer to the interpolation class
    autoPtr<interpolationCellPoint<vector>> Einterp_;

public:

    //- Runtime type information
    TypeName("VolumeWeighting");


    //- Constructors

        //- Construct from components
        VolumeWeightingFW
        (
            const word& fieldName,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        VolumeWeightingFW(const VolumeWeightingFW& am);

        //- Construct and return a clone
        virtual autoPtr<FieldWeigthing> clone() const
        {
            return autoPtr<FieldWeigthing>
            (
                new VolumeWeightingFW(*this)
            );
        }


    //- Destructor
    virtual ~VolumeWeightingFW();


    //- Member Functions

        //- Return electric field vector for the coordinates
        vector getFieldVector
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- Update the interpolation
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
