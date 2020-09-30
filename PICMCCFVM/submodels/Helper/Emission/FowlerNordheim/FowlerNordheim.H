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
    Foam::FowlerNordheim

Description
    Fowler Nordheim emission model

    R. H. Fowler & L. Nordheim
    Electron Emission in Intense Electric Field
    Proceedings of the Royal Society of London. Series A, Containing Papers of aMathematical and Physical Character, Vol. 119, No. 781 (May 1, 1928), pp. 173-181

SourceFiles
    FowlerNordheim.C

\*---------------------------------------------------------------------------*/

#ifndef FowlerNordheim_H
#define FowlerNordheim_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FowlerNordheim Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FowlerNordheim: public EmissionModel<CloudType>
{
    //- Parameter of the models
    scalar A_;
    scalar B_;
    scalar t_;

    //- Work function
    scalar fnW_;

    //- nEquivalentParticles of the emitted particles
    scalar nParticle_;

    //- Emission remainder
    scalar pRemainder_;
public:

    //- Runtime type information
    TypeName("FowlerNordheim");


    FowlerNordheim
    (
        const dictionary& dict,
        CloudType& cloud
    );

    virtual autoPtr<EmissionModel<CloudType>> clone() const
        {
            return autoPtr<EmissionModel<CloudType>>
            (
                new FowlerNordheim<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~FowlerNordheim();

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
    #include "FowlerNordheim.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
