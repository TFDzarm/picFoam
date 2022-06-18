/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2022 picFoam
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
    Foam::NoPush

Description
    No update of particle velocities.

\*---------------------------------------------------------------------------*/

#ifndef NoPush_H
#define NoPush_H

#include "ParticlePusher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class NoPush Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPush
:
    public ParticlePusher<CloudType>
{
//Private data

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoPush
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~NoPush();


    // Member Functions

        //- Do no parcels velocity update
        virtual void updateVelocity(typename CloudType::parcelType& p, const scalar trackTime) {}

        //- Return the unchanged velocity
        virtual vector correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const { return p.U(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoPush.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
