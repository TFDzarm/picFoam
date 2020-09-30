/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::NoElectronNeutralCollision

Description
    No collison.

\*---------------------------------------------------------------------------*/

#ifndef NoBinaryElectronCollision_H
#define NoBinaryElectronCollision_H

#include "ElectronNeutralCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoBinaryElectronCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoElectronNeutralCollision
:
    public ElectronNeutralCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoElectronNeutralCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~NoElectronNeutralCollision();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Do nothing
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, typename CloudType::parcelType& pN) {}

        //- Do nothing
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, vector& Un, label idN) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoElectronNeutralCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
