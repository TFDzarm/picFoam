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
    Foam::FloatingPotential

Description
    Floating potential boundary condition.

    Based on:
    J. P. Verboncoeur et al.
    Simultaneous potential and circuit solution for bounded plasma particle simulation codes
    https://ui.adsabs.harvard.edu/abs/1990spcs.rept.....V/abstract

\*---------------------------------------------------------------------------*/

#ifndef FloatingPotential_H
#define FloatingPotential_H

#include "BoundaryModel.H"
#include "fvMesh.H"
#include "circuitBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class FloatingPotential Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FloatingPotential
:
    public BoundaryModel<CloudType>
{
// Private data

    //- Charge gathered during current timestep
    scalar Qconv_;

    //- Surface charge
    scalar sigma_;

    //- Associated patchId
    label patchId_;

    //- Associated area
    scalar patchArea_;

    //- Boundary condition
    circuitBoundaryFvPatchField* boundaryCondition_;


public:

    //- Runtime type information
    TypeName("FloatingPotential");


    // Constructors

        //- Construct from dictionary
        FloatingPotential
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );

        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new FloatingPotential<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~FloatingPotential();


    // Member Functions

        //- Introduce particles
        virtual void injection();

        //- Called before and after Maxwell solver
        virtual void preUpdate_Boundary();
        virtual void postUpdate_Boundary();

        //- Patch boundary condition (e.g. reflect, delete, ...), return true if interaction occurred
        virtual bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        //- Called if particle was ejected on the patch (also called by BoundaryEvent)
        virtual void particleEjection(typename CloudType::parcelType& p, label patchId);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FloatingPotential.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
