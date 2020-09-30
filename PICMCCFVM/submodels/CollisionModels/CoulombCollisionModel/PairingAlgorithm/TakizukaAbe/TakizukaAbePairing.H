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
    Foam::TakizukaAbePairing

Description
    Pairing Algorithm: T. Takizuka & H. Abe
    A binary collision model for plasma simulation with a particle code
    https://doi.org/10.1016/0021-9991(77)90099-7

\*---------------------------------------------------------------------------*/

#ifndef TakizukaAbePairing_H
#define TakizukaAbePairing_H

#include "PairingAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class TakizukaAbePairing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TakizukaAbePairing
:
    public PairingAlgorithm<CloudType>
{
public:

    //- Runtime type information
    TypeName("TakizukaAbe");


    // Constructors

        //- Construct from dictionary
        TakizukaAbePairing
        (
            const dictionary& dict,
            CloudType& cloud,
            CoulombCollisionModel<CloudType>& collisionModel
        );


    //- Destructor
    virtual ~TakizukaAbePairing();


    // Member Functions

        //- Select particle pairs and collide
        virtual void pairANDcollide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TakizukaAbePairing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
