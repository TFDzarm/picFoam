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
    Foam::CountInsertions

Description
    Count insertions of particles into the domain.

SourceFiles
    CountInsertions.C

\*---------------------------------------------------------------------------*/

#ifndef CountInsertions_H
#define CountInsertions_H

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CountInsertions Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CountInsertions :
        public BoundaryEvent<CloudType>
{
    //Private data

        //- Number of particles inserted in the current timestep
        List<Field<scalar>> countParticle_;

        //- Number of parcels inserted in the current timestep
        List<Field<label>> countParcel_;

        //- Running average of inserted particles
        List<Field<scalar>> runningAverage_;

        //- Time when averaging began
        scalar averageStart_;
public:

    //- Runtime type information
    TypeName("CountInsertions");


    // Constructors

    CountInsertions(const CountInsertions<CloudType>& im);

    //- Construct from dictionary
        CountInsertions
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~CountInsertions();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new CountInsertions<CloudType>(*this)
        );
    }

    // Member Functions

        //- Do nothing with this
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td) {}

        //- Insertion
        virtual void ejectionEvent(typename CloudType::parcelType& p, label patchId);

        //- Print statistics
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CountInsertions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
