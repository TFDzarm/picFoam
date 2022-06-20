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
    Foam::UnidirectionalVelocityInfo

Description
    Calculates and prints info on particle velocities squared (<v^2>) based on given direction.

SourceFiles
    UnidirectionalVelocityInfo.C

\*---------------------------------------------------------------------------*/

#ifndef UnidirectionalVelocityInfo_H
#define UnidirectionalVelocityInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class UnidirectionalVelocityInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class UnidirectionalVelocityInfo :
        public DiagnosticInfo<CloudType>
{
    //- Direction mask to constrain velocity directions e.g. (0 1 1)
    vector directionMask_;

    //- Constrained velocity & velocity
    List<vector> velocitySqr_;

    //- Number of particles for each species
    List<scalar> nParticleTypes_;

public:

    //- Runtime type information
    TypeName("UnidirectionalVelocitySqr");


    // Constructors

    UnidirectionalVelocityInfo(const UnidirectionalVelocityInfo<CloudType>& im);

    //- Construct from dictionary
        UnidirectionalVelocityInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~UnidirectionalVelocityInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new UnidirectionalVelocityInfo<CloudType>(*this)
        );
    }
    // Member Functions

        //Called by DiagnosticsList so we loop only once through all particles
        virtual void gatherDiagnostic(const typename CloudType::parcelType& p);

        //- Print info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnidirectionalVelocityInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
