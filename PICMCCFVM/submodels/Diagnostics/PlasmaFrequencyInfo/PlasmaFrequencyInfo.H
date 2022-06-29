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
    Foam::PlasmaFrequencyInfo

Description
    Calculates and prints the plasma frequency of the electron species.

SourceFiles
    PlasmaFrequencyInfo.C

\*---------------------------------------------------------------------------*/

#ifndef PlasmaFrequencyInfo_H
#define PlasmaFrequencyInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PlasmaFrequencyInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PlasmaFrequencyInfo :
        public DiagnosticInfo<CloudType>
{
    //- Number of electron particles
    scalar nParticle_;

    //- Volume of the entire system
    scalar totalVolume_;

    //- Warn about a too large time step
    Switch warnTimeStep_;
public:

    //- Runtime type information
    TypeName("PlasmaFrequency");


    // Constructors

    PlasmaFrequencyInfo(const PlasmaFrequencyInfo<CloudType>& im);

    //- Construct from dictionary
        PlasmaFrequencyInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PlasmaFrequencyInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new PlasmaFrequencyInfo<CloudType>(*this)
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
    #include "PlasmaFrequencyInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
