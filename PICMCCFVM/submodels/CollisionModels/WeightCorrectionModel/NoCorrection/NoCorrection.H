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
    Foam::NoWeightCorrectionModel

Description
    Do not correct the velocities.

SourceFiles
    NoCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef NoWeightCorrectionModel_H
#define NoWeightCorrectionModel_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "WeightCorrectionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class WallInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoWeightCorrectionModel
:
        public WeightCorrectionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoWeightCorrectionModel
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~NoWeightCorrectionModel();


    // Member Functions

        //- Correct the velocity
        virtual void correctVelocity(
                typename CloudType::parcelType* parcelP,
                typename CloudType::parcelType* parcelQ,
                vector preUp,
                vector preUq
                );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "NoCorrection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
