/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::ZARMInOutflow

Description
    Modification to FreeStream, reflect particles at the inlet.

\*---------------------------------------------------------------------------*/

#ifndef ZARMInOutflow_H
#define ZARMInOutflow_H

#include "BoundaryModel.H"
#include "volFields.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class ZARMInOutflow Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ZARMInOutflow
:
    public BoundaryModel<CloudType>
{
    // Private data

        //- Patch type
        enum patchType
        {
            pTNoInteraction,
            pTInletOutlet,
            pTInlet
        };

        //- The indices of patches to introduce molecules across
        labelList patches_;

        //- The molecule types to be introduced
        List<List<label>> moleculeTypeIds_;

        //- The number density of the species in the inflow for every patch
        List<Field<scalar>> numberDensities_;

        //- A List of Lists of Fields specifying carry-over of mass flux from
        // one timestep to the next
        // + Outer List - one inner List for each patch
        // + Inner List - one Field for every species to be introduced
        // + Each field entry corresponding to a face to be injected across
        //   with a particular species
        List<List<Field<scalar>>> particleFluxAccumulators_;

        //- List of patch types for all patches
        List<patchType> interactionList_;

        //- Boundary temperatures for each species separately
        List<Field<scalar>> boundaryT_;

        //- Boundary velocity for each species separately
        List<Field<vector>> boundaryU_;
public:

    //- Runtime type information
    TypeName("ZARMInOutflow");


    // Constructors

        //- Construct from dictionary
        ZARMInOutflow
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new ZARMInOutflow<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~ZARMInOutflow();


    // Member Functions

        //- Introduce particles
        virtual void injection();

        bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        using BoundaryModel<CloudType>::particleEjection;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ZARMInOutflow.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
