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
    Foam::NoCoulombCollision

Description
    No coulomb collison model.

\*---------------------------------------------------------------------------*/

#ifndef NoCoulombCollision_H
#define NoCoulombCollision_H

#include "CoulombCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoCoulombCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoCoulombCollision
:
    public CoulombCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoCoulombCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~NoCoulombCollision();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Do nothing...
        virtual void collide
        (
            typename CloudType::parcelType* pP,
            typename CloudType::parcelType* pQ,
            scalar debyeLength,
            scalar nP,
            scalar nQ,
            scalar nPQ
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoCoulombCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
