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
    Foam::EquipartitionInitialization

Description
    Equipartitioned initialization.

\*---------------------------------------------------------------------------*/

#ifndef EquipartitionInitialization_H
#define EquipartitionInitialization_H

#include "InitializationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class EquipartitionInitialization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EquipartitionInitialization
:
    public InitializationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("Equipartition");


    // Constructors

        //- Construct from dictionary
        EquipartitionInitialization
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~EquipartitionInitialization();


    // Member Functions

        //- Ininitalize ...
        virtual void initialiseParticles();

        //- Add particle called in initialiseParticles()
        void addParticles(label celli, const vector& velocity);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EquipartitionInitialization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

