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
    Foam::HigueraCaryPusher

Description
    Particle pusher by A.V. Higuera & J.R. Cary

    Structure-preserving second-order integration of relativistic charged particle
    trajectories in electromagnetic fields
    https://doi.org/10.1063/1.4979989

\*---------------------------------------------------------------------------*/

#ifndef HigueraCary_H
#define HigueraCary_H

#include "ParticlePusher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class HigueraCaryPusher Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HigueraCaryPusher
:
    public ParticlePusher<CloudType>
{
    //- If there is no magnetic field, skip rotation procedure.
    bool calcMagneticRotation_;
public:

    //- Runtime type information
    TypeName("HigueraCary");


    // Constructors

        //- Construct from dictionary
        HigueraCaryPusher
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~HigueraCaryPusher();


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
    #include "HigueraCary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
