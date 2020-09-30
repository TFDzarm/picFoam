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
    Foam::ElectroStaticSolver

Description
    Electrostatic field solver. Solves Poisson's equation.
\*---------------------------------------------------------------------------*/

#ifndef ElectroStatic_H
#define ElectroStatic_H

#include "MaxwellSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class ElectroStaticSolver Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ElectroStaticSolver
:
    public MaxwellSolver<CloudType>
{
public:

    //- Runtime type information
    TypeName("ElectroStatic");


    // Constructors

        //- Construct from dictionary
        ElectroStaticSolver
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~ElectroStaticSolver();


    // Member Functions

        //- Solve Poisson equation
        virtual void solveFields();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ElectroStaticSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
