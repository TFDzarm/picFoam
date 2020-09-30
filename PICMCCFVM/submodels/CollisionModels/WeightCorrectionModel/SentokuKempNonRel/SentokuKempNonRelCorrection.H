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
    Foam::SentokuKempNonRelCorrection

Description

    Correction model based on

    Y. Sentoku & A.J. Kemp

    Numerical methods for particle simulations at extreme densities and temperatures: Weighted particles, relativistic collisions and reduced currents
    https://doi.org/10.1016/j.jcp.2008.03.043

    Here we perform calculation nonrelativistic!

    WARNING:    Notice an error occurs when the heavier species has a lower weight => the lighter species is corrected (has more particles in on parcel).
                The other way arround it works fine! Implementation is exactly as the paper states and is the same as in EPOCH. Model is incomplete/ not correct?
                ToDo: Investigate

    Code structure inspired by:
    EPOCH
    Code:  https://cfsa-pmw.warwick.ac.uk/users/sign_in
    Paper: https://doi.org/10.1088%2F0741-3335%2F57%2F11%2F113001

SourceFiles
    SentokuKempNonRelCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef SentokuKempNonRelCorrection_H
#define SentokuKempNonRelCorrection_H

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
class SentokuKempNonRelCorrection
:
        public WeightCorrectionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("SentokuKempNR");

    // Constructors

        //- Construct from components
        SentokuKempNonRelCorrection
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~SentokuKempNonRelCorrection();


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
    #include "SentokuKempNonRelCorrection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
