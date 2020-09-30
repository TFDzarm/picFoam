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
    Foam::NanbuCollision

Description
    Collision model by: K. Nanbu
    Probability Theory of Electron–Molecule, Ion–Molecule, Molecule–Molecule, and Coulomb Collisions for Particle Modeling of Materials Processing Plasmas and Gases
    https://doi.org/10.1109/27.887765

\*---------------------------------------------------------------------------*/

#ifndef NanbuCollision_H
#define NanbuCollision_H

#include "CoulombCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NanbuCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NanbuCollision
:
    public CoulombCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("Nanbu");


    // Constructors

        //- Construct from dictionary
        NanbuCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~NanbuCollision();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Apply collision
        virtual void collide
        (
            typename CloudType::parcelType* parcelP,
            typename CloudType::parcelType* parcelQ,
            scalar debyeLength,
            scalar nP,
            scalar nQ,
            scalar nPQ
        );

        //- Calculate scattering angle
        scalar calculate_cosChi(scalar s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NanbuCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
