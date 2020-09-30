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
    Foam::TemperatureInfo

Description
    Calculates and prints info on temperature.

SourceFiles
    TemperatureInfo.C

\*---------------------------------------------------------------------------*/

#ifndef TemperatureInfo_H
#define TemperatureInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TemperatureInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TemperatureInfo :
        public DiagnosticInfo<CloudType>
{
    //- Velocity <v^2>
    List<scalar> vSqr_;

    //- Velocity <v>
    List<vector> vDrift_;

    //- Number of particles for each species
    List<scalar> nParticleTypes_;

    //- Include drift velocity in calculation?
    Switch accountForDrift_;

    //- Total number of particles
    scalar nMoles_;
public:

    //- Runtime type information
    TypeName("Temperature");


    // Constructors

    TemperatureInfo(const TemperatureInfo<CloudType>& im);

    //- Construct from dictionary
        TemperatureInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~TemperatureInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new TemperatureInfo<CloudType>(*this)
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
    #include "TemperatureInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
