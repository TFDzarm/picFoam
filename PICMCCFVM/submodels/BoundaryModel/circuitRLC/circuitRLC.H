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
    Foam::circuitRLC

Description
    RLC circuit model. Translated from the original finite difference model to finite volume model

    Based on:
    J. P. Verboncoeur et al.
    Simultaneous potential and circuit solution for bounded plasma particle simulation codes
    https://ui.adsabs.harvard.edu/abs/1990spcs.rept.....V/abstract

\*---------------------------------------------------------------------------*/

#ifndef circuitRLC2_H
#define circuitRLC2_H

#include "BoundaryModel.H"
#include "fvMesh.H"
#include "circuitBoundaryFvPatchField.H"
#include "RLCBoundaryFvPatchField.H"
#include "EmissionModelList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class circuitRLC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class circuitRLC
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
    scalar sigma1_;

    //- Charge gathered during current timestep
    scalar chargeAccumulator_;


    //- Cathode patch id
    label cathodePatch_;

    //- Anode patch id
    label anodePatch_;

    //- Cathode patch area
    scalar cathodeArea_;

    //RLC boundaty on the cathode
    RLCBoundaryFvPatchField* boundaryCondition_;

    //- Provided potential
    scalar V_;

    //- Used for sinusoid voltage source
    scalar omega_;
    scalar phase_;

    //- Charges in the circuit
    scalar Qex_; //current timestep
    scalar Qex1_; //timestep-1
    scalar Qex2_; //-2
    scalar Qex3_; //-3

    //- Resistance
    scalar R_;

    //- Capacity
    scalar C_;

    //- Inductance
    scalar L_;

    //- Integration parameter
    scalar a0_, a1_, a2_, a3_, a4_;
    scalar k_;

    //- Circuit current
    scalar Iex_;

    //- Permittivity (can be changed for testing purposes)
    scalar epsilon_;

    //Emission models
    EmissionModelList<CloudType> emissionList_;
public:

    //- Runtime type information
    TypeName("circuitRLC");


    // Constructors

        //- Construct from dictionary
        circuitRLC
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );

        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new circuitRLC<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~circuitRLC();


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

        //- Return external current
        inline scalar current() const;

        //- Return average of the boundary potential
        inline scalar voltage() const;

        inline scalar source() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "circuitRLCI.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "circuitRLC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
