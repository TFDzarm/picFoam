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
    Foam::PerezRelativisticCollision

Description
     Collision model from: F. Pérez et al.

     Improved modeling of relativistic collisions and collisional ionization in particle-in-cell codes
     https://doi.org/10.1063/1.4742167

\*---------------------------------------------------------------------------*/

#ifndef PerezRelativisticCollision_H
#define PerezRelativisticCollision_H

#include "CoulombCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class PerezRelativisticCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PerezRelativisticCollision
:
    public CoulombCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("Perez");


    // Constructors

        //- Construct from dictionary
        PerezRelativisticCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PerezRelativisticCollision();


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
    #include "PerezRelativisticCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
