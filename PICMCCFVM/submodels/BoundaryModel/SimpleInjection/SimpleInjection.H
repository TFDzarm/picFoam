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
    Foam::SimpleInjection

Description
    SimpleInjection model. Inject particles with given injection rate.
    Boundaries can be set to be reflecting, floating or simiply deleting.

\*---------------------------------------------------------------------------*/

#ifndef SimpleInjection_H
#define SimpleInjection_H

#include "BoundaryModel.H"
#include "fvMesh.H"

#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "faceList.H"
#include "ParticleEmitter.H"
#include "circuitBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class SimpleInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SimpleInjection
:
    public BoundaryModel<CloudType>
{
// Private data

    //- Patch types
    enum patchType
    {
        ptDeletionPatch,
        ptReflection,
        ptFloating
    };

    //- Emitter helper class
    ParticleEmitter<CloudType> emitter_;

    //- PatchId of the injector patch
    label injectorPatch_;

    //- The types to be introduced
    List<word> injectionSpeciesList_;
    List<label> injectionSpeciesTypeIds_;

    //- List of all patches (excluding processor patches) and their type
    List<patchType> interactionList_;

    //Injecton Frequency [1/s]
    List<scalar> iF_;
    List<scalar> particleRemainder_;

    //-Temperatur [eV]
    List<scalar> eV_;

    //- If we reflect should we reflect with temperature according to eV_?
    List<bool> reflectWithInitVel_;

    //- Use a floating boundary condition?
    bool isFloatingBC_;

    //- Charge gathered during current timestep
    scalar Qconv_;

    //- Patch area
    scalar patchArea_;

    //- Boundary condition
    circuitBoundaryFvPatchField* boundaryCondition_;

    //- Time at which no particles should be created
    scalar injectionEndTime_;

public:

    //- Runtime type information
    TypeName("SimpleInjection");


    // Constructors

        //- Construct from dictionary
        SimpleInjection
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );

        virtual autoPtr<BoundaryModel<CloudType>> clone() const
        {
            return autoPtr<BoundaryModel<CloudType>>
            (
                new SimpleInjection<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~SimpleInjection();


    // Member Functions

        //- Introduce particles
        virtual void injection();

        //- Called before Maxwell solver
        virtual void preUpdate_Boundary();

        //- Patch boundary condition (e.g. reflect, delete, ...), return true if interaction occurred
        bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        //- Called if particle was ejected on the patch (also called by BoundaryEvent)
        void particleEjection(typename CloudType::parcelType& p, label patchId);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SimpleInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
