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
    Foam::RelativisticElectronNeutralCollision

Description
    Relativistic collisions of electron and neutrals.
    Model is based on:

    F. Pérez et al.
    Improved modeling of relativistic collisions and collisional ionization in particle-in-cell codes
    https://doi.org/10.1063/1.4742167

    and uses the scattering angle calulation of

    K. Nanbu
    Probability Theory of Electron–Molecule, Ion–Molecule, Molecule–Molecule, and Coulomb Collisions for Particle Modeling of Materials Processing Plasmas and Gases
    https://doi.org/10.1109/27.887765

\*---------------------------------------------------------------------------*/

#ifndef RelativisticElectronNeutralCollision_H
#define RelativisticElectronNeutralCollision_H

#include "ElectronNeutralCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class RelativisticElectronNeutralCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RelativisticElectronNeutralCollision
:
    public ElectronNeutralCollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("RelativisticElectronNeutral");


    // Constructors

        //- Construct from dictionary
        RelativisticElectronNeutralCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~RelativisticElectronNeutralCollision();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Update the velocities
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, typename CloudType::parcelType& pN);

        //- Update the velocities for collisions with the background
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, vector& Un, label idN);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RelativisticElectronNeutralCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
