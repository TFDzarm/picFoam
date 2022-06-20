/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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
    Foam::VelocityInfo

Description
    Calculates and prints info on the drift velocity.

SourceFiles
    VelocityInfo.C

\*---------------------------------------------------------------------------*/

#ifndef VelocityInfo_H
#define VelocityInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class VelocityInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VelocityInfo :
        public DiagnosticInfo<CloudType>
{
    //- Velocity <v>
    List<vector> vDrift_;

    //- Number of particles for each species
    List<scalar> nParticleTypes_;

public:

    //- Runtime type information
    TypeName("Velocity");


    // Constructors

    VelocityInfo(const VelocityInfo<CloudType>& im);

    //- Construct from dictionary
        VelocityInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~VelocityInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new VelocityInfo<CloudType>(*this)
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
    #include "VelocityInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
