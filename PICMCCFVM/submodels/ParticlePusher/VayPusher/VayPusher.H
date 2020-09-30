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
    Foam::VayPusher

Description
    Particle pusher by J.-L. Vay

    Simulation of beams or plasmas crossing at relativistic velocity
    https://doi.org/10.1063/1.2837054

\*---------------------------------------------------------------------------*/

#ifndef VayPusher_H
#define VayPusher_H

#include "ParticlePusher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class VayPusher Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VayPusher
:
    public ParticlePusher<CloudType>
{
public:

    //- Runtime type information
    TypeName("Vay");


    // Constructors

        //- Construct from dictionary
        VayPusher
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~VayPusher();


    // Member Functions

        //- Update the parcels velocity
        virtual void updateVelocity(typename CloudType::parcelType& p, const scalar trackTime);

        //- Return the velocity change
        virtual vector correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VayPusher.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
