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
    Foam::VariableHardSphere

Description
    Variable hard sphere model.

\*---------------------------------------------------------------------------*/

#ifndef VariableHardSphere_H
#define VariableHardSphere_H

#include "TotalCrossSectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class VariableHardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VariableHardSphere
:
    public TotalCrossSectionModel<CloudType>
{
    // Private data

        //- Reference temperature
        const scalar Tref_;

public:

    //- Runtime type information
    TypeName("VariableHardSphere");


    // Constructors

        //- Construct from dictionary
        VariableHardSphere
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~VariableHardSphere();


    // Member Functions

        //- Total cross section
        virtual scalar sigmaTcR
        (
            const typename CloudType::parcelType& pP,
            const typename CloudType::parcelType& pQ
        ) const;

        //- Total cross section (background)
        virtual scalar sigmaTcR
        (
            const typename CloudType::parcelType& pP,
            const vector& Uq,
            label idQ
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VariableHardSphere.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

