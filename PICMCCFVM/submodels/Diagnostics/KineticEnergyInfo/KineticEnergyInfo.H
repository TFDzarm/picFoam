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
    Foam::KineticEnergyInfo

Description
    Prints info on kinetic energy for all species.

SourceFiles
    KineticEnergyInfo.C

\*---------------------------------------------------------------------------*/

#ifndef KineticEnergyInfo_H
#define KineticEnergyInfo_H

#include "DiagnosticInfo.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class KineticEnergyInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KineticEnergyInfo :
        public DiagnosticInfo<CloudType>
{

    //- KineticEnergy of types
    List<scalar> kineticEnergyofTypes_;

    //- Calculate relativisticly?
    Switch relativistic_;
public:

    //- Runtime type information
    TypeName("KineticEnergy");


    // Constructors

    KineticEnergyInfo(const KineticEnergyInfo<CloudType>& im);

    //- Construct from dictionary
        KineticEnergyInfo
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~KineticEnergyInfo();


    virtual autoPtr<DiagnosticInfo<CloudType>> clone() const
    {
        return autoPtr<DiagnosticInfo<CloudType>>
        (
            new KineticEnergyInfo<CloudType>(*this)
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
    #include "KineticEnergyInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
