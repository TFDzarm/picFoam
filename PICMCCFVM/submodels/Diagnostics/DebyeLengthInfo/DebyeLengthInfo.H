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
    Foam::DebyeLengthInfo

Description
    Calculates and prints the average Debye length.

SourceFiles
    DebyeLengthInfo.C

\*---------------------------------------------------------------------------*/

#ifndef DebyeLengthInfo_H
#define DebyeLengthInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DebyeLengthInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DebyeLengthInfo :
        public DiagnosticInfo<CloudType>
{
    //- Velocity <v^2>
    List<scalar> vSqr_;

    //- Velocity <v>
    List<vector> vDrift_;

    //- Number of particles for each species
    List<scalar> nParticleTypes_;

    //- Charge <q>
    List<scalar> sCharge_;

    //- Volume of the entire system
    scalar totalVolume_;

    //- Include drift velocity in calculation?
    Switch accountForDrift_;

    //- Global average cell length scale
    scalar gCellLengthScale_;

    //- Warn about a too large cell size
    Switch warnCellSize_;
public:

    //- Runtime type information
    TypeName("DebyeLength");


    // Constructors

    DebyeLengthInfo(const DebyeLengthInfo<CloudType>& im);

    //- Construct from dictionary
        DebyeLengthInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~DebyeLengthInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new DebyeLengthInfo<CloudType>(*this)
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
    #include "DebyeLengthInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
