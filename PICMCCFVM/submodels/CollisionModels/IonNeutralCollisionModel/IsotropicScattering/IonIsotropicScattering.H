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
    Foam::IsotropicScattering

Description
    Simple isotropic scattering model. Used to update the velocities of particles colliding.

\*---------------------------------------------------------------------------*/

#ifndef IonIsotropicScattering_H
#define IonIsotropicScattering_H

#include "IonNeutralCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class IsotropicScattering Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IonIsotropicScattering
:
    public IonNeutralCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("IsotropicScattering");


    // Constructors

        //- Construct from dictionary
        IonIsotropicScattering
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~IonIsotropicScattering();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;


        //- Update the velocity
        virtual void updateVelocity
        (
            typename CloudType::parcelType& pP,
            typename CloudType::parcelType& pQ
        );

        //- Update the velocity (for background collisions)
        virtual void updateVelocity
        (
            typename CloudType::parcelType& pP,
            vector& Uq,
            label idQ
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IonIsotropicScattering.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
