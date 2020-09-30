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
    Foam::IdealCurrentSource

Description
    IdealCurrent boundary condition.

    Based on:
    J. P. Verboncoeur et al.
    Simultaneous potential and circuit solution for bounded plasma particle simulation codes
    https://ui.adsabs.harvard.edu/abs/1990spcs.rept.....V/abstract

\*---------------------------------------------------------------------------*/

#ifndef IdealCurrentSource_H
#define IdealCurrentSource_H

#include "BoundaryModel.H"
#include "fvMesh.H"
#include "circuitBoundaryFvPatchField.H"
#include "EmissionModelList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class FreeStream Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IdealCurrentSource
:
    public BoundaryModel<CloudType>
{
// Private data

    //- Patch types
    enum electrodeType
    {
        etNone,
        etCathode,
        etAnode
    };

    //- List of all patches (excluding processor patches) and their type
    List<electrodeType> interactionList_;


    //- Surface charge
    scalar sigma_;

    //- Charge gathered during current timestep
    scalar chargeAccumulator_;

    //- Cathode patch id
    label cathodePatch_;

    //- Anode patch id
    label anodePatch_;

    //- Anode patch area
    scalar anodeArea_;

    //- Cathode patch area
    scalar cathodeArea_;

    //Boundary condition
    circuitBoundaryFvPatchField* boundaryCondition_;

    //- Total current supplied to the cathode
    scalar I_;

    //- Total current density supplied to the cathode
    scalar J_;

    //- Emission models
    EmissionModelList<CloudType> emissionList_;
public:

    //- Runtime type information
    TypeName("IdealCurrentSource");


    // Constructors

        //- Construct from dictionary
        IdealCurrentSource
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );

        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new IdealCurrentSource<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~IdealCurrentSource();


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
    #include "IdealCurrentSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
