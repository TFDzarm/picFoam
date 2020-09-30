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
    Foam::NoIonization

Description
    No ionization in coulomb collisions.

SourceFiles
    NoIonization.C

\*---------------------------------------------------------------------------*/

#ifndef NoIonization_H
#define NoIonization_H

#include "IonizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NoIonization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoIonization :
        public IonizationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from dictionary
        NoIonization
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~NoIonization();


    // Member Functions

        //- Initialize
        virtual void initIonization(label pTypeId, label qTypeId, scalar nP, scalar nQ);

        //- Prepare density ratios before ionizations
        virtual void prepareNumberDensities(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ);

        //- Perpare ionization of the parcel
        virtual void prepareIonization(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ);

        //- Ionize the parcel
        virtual bool ionize();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoIonization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
