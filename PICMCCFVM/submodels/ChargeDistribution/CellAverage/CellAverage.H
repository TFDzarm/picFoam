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
    Foam::CellAverage

Description
    Add the particles charges to the cell centered field and divide by the cells volume.

\*---------------------------------------------------------------------------*/

#ifndef CellAverage_H
#define CellAverage_H

#include "ChargeDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CellAverage Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CellAverage :
        public ChargeDistribution<CloudType>
{
public:

    //- Runtime type information
    TypeName("CellAverage");


    // Constructors

    //- Construct from dictionary
        CellAverage
        (
            const word& fieldName,
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~CellAverage();


    // Member Functions

        //- Add particle charge to the field
        void add(const typename CloudType::parcelType& parcel);

        //- Update the distribution
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CellAverage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
