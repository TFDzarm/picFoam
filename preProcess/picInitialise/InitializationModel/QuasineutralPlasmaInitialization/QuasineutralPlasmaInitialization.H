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
    Foam::QuasineutralPlasmaInitialization

Description
    Initalize a quasineutral plasma

\*---------------------------------------------------------------------------*/

#ifndef QuasineutralPlasmaInitialization_H
#define QuasineutralPlasmaInitialization_H

#include "InitializationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class QuasineutralPlasmaInitialization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class QuasineutralPlasmaInitialization
:
    public InitializationModel<CloudType>
{
//Private data

    //- The ion species
    word ionSpecies_;

    //- TypeId of the ion species
    label ionTypeId_;

    //- Number of particles per cell so that the plasma is quasineutral
    List<label> cellCount_;

    //- Calc velocity from maxwellian distribution or mean kinetic energy
    bool useMaxwellian_;

public:

    //- Runtime type information
    TypeName("QuasineutralPlasma");


    // Constructors

        //- Construct from dictionary
        QuasineutralPlasmaInitialization
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~QuasineutralPlasmaInitialization();


    // Member Functions

        //- Initialize...
        virtual void initialiseParticles();

        //- Set by model do we read before or calculate after initialization
        virtual bool calculateNumberDensities()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QuasineutralPlasmaInitialization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

