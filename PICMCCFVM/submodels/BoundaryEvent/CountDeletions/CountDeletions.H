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
    Foam::CountDeletions

Description
    Count deletions of particles leaving domain.

SourceFiles
    CountDeletions.C

\*---------------------------------------------------------------------------*/

#ifndef CountDeletions_H
#define CountDeletions_H

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CountDeletions Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CountDeletions :
        public BoundaryEvent<CloudType>
{
    //Private data

        //- Number of particles delted in the current timestep
        List<Field<scalar>> countParticle_;

        //- Number of parcels delted in the current timestep
        List<Field<label>> countParcel_;

        //- Running average of deleted particles
        List<Field<scalar>> runningAverage_;

        //- Time when averaging began
        scalar averageStart_;
public:

    //- Runtime type information
    TypeName("CountDeletions");


    // Constructors

    CountDeletions(const CountDeletions<CloudType>& im);

    //- Construct from dictionary
        CountDeletions
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~CountDeletions();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new CountDeletions<CloudType>(*this)
        );
    }

    // Member Functions

        //- Check deletion
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        //- Print statistics
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CountDeletions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
