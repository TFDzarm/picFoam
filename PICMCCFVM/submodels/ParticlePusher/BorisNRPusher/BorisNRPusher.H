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
    Foam::BorisNRPusher

Description
    Nonrelativstic version of Boris pusher.

    J.P. Boris
    Relativistic plasma simulation — Optimization of a hybrid code
    Proceedings of 4th Conference on Numerical Simulation of Plasmas 1970
    pp. 3-67

\*---------------------------------------------------------------------------*/

#ifndef BorisNRPusher_H
#define BorisNRPusher_H

#include "ParticlePusher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class BorisNRPusher Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BorisNRPusher
:
    public ParticlePusher<CloudType>
{
//Private data

    //- If there is no magnetic field, skip rotation procedure.
    bool calcMagneticRotation_;
public:

    //- Runtime type information
    TypeName("BorisNR");


    // Constructors

        //- Construct from dictionary
        BorisNRPusher
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~BorisNRPusher();


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
    #include "BorisNRPusher.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
