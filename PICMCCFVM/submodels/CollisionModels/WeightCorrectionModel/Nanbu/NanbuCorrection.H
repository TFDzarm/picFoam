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
    Foam::NanbuWeightCorrectionModel

Description
    Correction for collision of particles with unequal weight.

    K.Nanbu & S.Yonemura
    Weighted Particles in Coulomb Collision Simulations Based on the Theory of a Cumulative Scattering Angle
    https://doi.org/10.1006/jcph.1998.6049

SourceFiles
    NanbuCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef NanbuWeightCorrectionModel_H
#define NanbuWeightCorrectionModel_H

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
class NanbuWeightCorrectionModel
:
        public WeightCorrectionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("Nanbu");

    // Constructors

        //- Construct from components
        NanbuWeightCorrectionModel
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~NanbuWeightCorrectionModel();


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
    #include "NanbuCorrection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
