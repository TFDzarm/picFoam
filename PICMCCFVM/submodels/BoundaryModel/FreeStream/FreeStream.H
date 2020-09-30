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
    Foam::FreeStream

Description
    Inserting new particles across the faces of a all patched of type
    "patch" for a free stream.  Uniform values number density, temperature
    and velocity sourced face-by-face from the boundaryT and boundaryU fields
    of the cloud.

    Modified OpenFOAM DSMC implementaion to fit in this framework.
\*---------------------------------------------------------------------------*/

#ifndef FreeStream_H
#define FreeStream_H

#include "BoundaryModel.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class FreeStream Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FreeStream
:
    public BoundaryModel<CloudType>
{
    // Private data

        //- The indices of patches to introduce molecules across
        labelList patches_;

        //- The molecule types to be introduced
        List<label> moleculeTypeIds_;

        //- The number density of the species in the inflow
        Field<scalar> numberDensities_;

        //- A List of Lists of Fields specifying carry-over of mass flux from
        // one timestep to the next
        // + Outer List - one inner List for each patch
        // + Inner List - one Field for every species to be introduced
        // + Each field entry corresponding to a face to be injected across
        //   with a particular species
        List<List<Field<scalar>>> particleFluxAccumulators_;


public:

    //- Runtime type information
    TypeName("FreeStream");


    // Constructors

        //- Construct from dictionary
        FreeStream
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new FreeStream<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~FreeStream();


    // Member Functions

        //- Introduce particles
        virtual void injection();

        using BoundaryModel<CloudType>::particleEjection;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FreeStream.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
