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
    Foam::ParticleEmitter

Description
    Particle emitter helper class

SourceFiles
    ParticleEmitter.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleEmitter_H
#define ParticleEmitter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "faceList.H"
#include "polyPatch.H"
#include "triFaceList.H"

namespace Foam
{

    template<class CloudType>
    class ParticleEmitter
    {
    // Private data

        //- Reference to cloud
        CloudType& cloud_;

        //- Patch id we emit from
        label patchId_;

        //- Patch area - total across all processors
        scalar patchArea_;

        //- Patch face normal directions
        vectorList patchNormals_;

        //- List of cell labels corresponding to injector positions
        labelList cellOwners_;

        //- Decomposed patch faces as a list of triangles
        triFaceList triFaces_;

        //- Addressing from per triangle to patch face
        labelList triToFace_;

        //- Cumulative triangle area per triangle face
        scalarList triCumulativeMagSf_;

        //- Cumulative area fractions per processor
        scalarList sumTriMagSf_;

        //- Function pointer to velocity model in use
        vector (ParticleEmitter::*getVelocity_)(const vector&, const vector&,
                                          const vector&, scalar, label);
    public:

        //Access

            //- Return const reference to cloud
            inline const CloudType& cloud() const;

            //- Return non-const reference to cloud
            inline CloudType& cloud();

            //- Return patch normal vector list
            inline const vectorList& patchNormals() const;

            //- Return patchId we emit from
            inline const label& patchId() const;

            //- Return patch area
            inline const scalar& patchArea() const;
    public:
        //Public data

            //- Available velocity models
            enum VelocityModel
            {
                vmUnknown = -1,
                vmBirdMaxwellianFlux,
                vmMostProbableSpeed,
                vmHalfMaxwellianFlux,
                vmMaxwellianFlux,
            };
    private:
        // Private data

            //-Choosen velocity model
            VelocityModel velocityModel_;
    public:
        //Constructor
            ParticleEmitter(CloudType& cloud);

            //- Set this in construction of model that uses this emiter!
            void initilizeParticleEmitter(label patchId, VelocityModel model);

        // Member functions

            //- Emit particle, choose random processor
            typename CloudType::parcelType* emitParticle(scalar temperature, label typeId, scalar nParticle = -1);

            //- Emit explicit on the processor this is called from
            typename CloudType::parcelType* emitParticleExplicit(scalar temperature, label typeId, scalar stepFraction, scalar nParticle = -1);

            //- Sample velocity from choosen model
            vector getVelocity(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId);
    private:
            //- Velocity models

            //-> vmBirdMaxwellianFlux
            vector getVelocity0(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId);

            //-> vmMostProbableSpeed
            vector getVelocity1(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId);

            //-> vmHalfMaxwellianFlux
            vector getVelocity2(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId);

            //-> vmMaxwellianFlux
            vector getVelocity3(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId);
    };

} // End namespace Foam

#include "ParticleEmitterI.H"

#ifdef NoRepository
    #include "ParticleEmitter.C"
#endif

#endif
