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
    Foam::VolumeWeighting

Description
    Use the barycentric coordinates of the particle to weight the charge to the cell center and vertices.
    Communicate vertex data and interpolate, data back to the cell center using the vertex cell-center distance. Divide by the cell volume.

\*---------------------------------------------------------------------------*/

#ifndef VolumeWeighting_H
#define VolumeWeighting_H

#include "ChargeDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class WallInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VolumeWeighting :
        public ChargeDistribution<CloudType>
{
//Private data

    //- Charge interpolated to the cell vertices
    Field<scalar> vertexData_;

    //- Charge interpolated to the cell center
    Field<scalar> ccData_;

    //- Weight for point to cellcenter interpolation
    scalarListList weights_;

public:

    //- Runtime type information
    TypeName("VolumeWeighting");


    // Constructors

    //- Construct from dictionary
        VolumeWeighting
        (
            const word& fieldName,
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~VolumeWeighting();


    // Member Functions

        //- Add particle charge to the field
        void add(const typename CloudType::parcelType& parcel);

        //- Update the distribution
        void update();

        //- Reset the charge density field
        void reset();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VolumeWeighting.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
