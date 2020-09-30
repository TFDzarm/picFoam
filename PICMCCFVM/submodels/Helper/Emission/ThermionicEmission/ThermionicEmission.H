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
    Foam::ThermionicEmission

Description
    Thermionic emission with a fixed temperature according to Richardson's law.

    C.R.Crowell
    The Richardson constant for thermionic emission in Schottky barrier diodes
    https://doi.org/10.1016%2F0038-1101%2865%2990116-4

    ToDo: Including the Schottky effect

    M.E. Kiziroglou et al.
    Thermionic field emission at electrodeposited Ni–Si Schottky barriers
    https://doi.org/10.1016%2Fj.sse.2008.03.002

SourceFiles
    ThermionicEmission.C

\*---------------------------------------------------------------------------*/

#ifndef ThermionicEmission_H
#define ThermionicEmission_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ThermionicEmission Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermionicEmission: public EmissionModel<CloudType>
{
    //- Model parametr
    scalar lambdaR_;
    scalar A0_;

    //- Fixed temperature
    scalar Tbc_;

    //- Work function
    scalar fnW_;

    //- Calculated current density
    scalar j_;

    //- nEquivalentParticles of the emitted particles
    scalar nParticle_;

    //- Emission remainder
    scalar pRemainder_;
public:

    //- Runtime type information
    TypeName("Thermionic");


    ThermionicEmission
    (
        const dictionary& dict,
        CloudType& cloud
    );

    virtual autoPtr<EmissionModel<CloudType>> clone() const
        {
            return autoPtr<EmissionModel<CloudType>>
            (
                new ThermionicEmission<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~ThermionicEmission();


    // Member Functions

        //- Emission of particles from boundary
        virtual void emission();

        //- Collisional emission e.g. sputter
        virtual void collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td){}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermionicEmission.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
