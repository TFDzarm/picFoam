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
    Foam::IsotropizationInitialization

Description
    Initalize perpendicular and parallel temperature components differently

\*---------------------------------------------------------------------------*/

#ifndef IsotropizationInitialization_H
#define IsotropizationInitialization_H

#include "InitializationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class IsotropizationInitialization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IsotropizationInitialization
:
    public InitializationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("Isotropization");


    // Constructors

        //- Construct from dictionary
        IsotropizationInitialization
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~IsotropizationInitialization();


    // Member Functions

        //- Set by model do we read before or calculate after initialization
        virtual void initialiseParticles();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IsotropizationInitialization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

