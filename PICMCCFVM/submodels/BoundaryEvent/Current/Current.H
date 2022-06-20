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
    Foam::Current

Description
    Calculate average current flowing over a patch

SourceFiles
    Current.C

\*---------------------------------------------------------------------------*/

#ifndef Current_H
#define Current_H

#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Current Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Current :
        public BoundaryEvent<CloudType>
{
    //Private data

        //-Charge flown past the patches in current timestep
        List<Field<scalar>> chargeCounter_;

        //- Running average of deleted particles
        List<Field<scalar>> runningAverage_;

        //- Time when averaging began
        scalar averageStart_;
public:

    //- Runtime type information
    TypeName("Current");


    // Constructors

    Current(const Current<CloudType>& im);

    //- Construct from dictionary
        Current
        (
            const dictionary& dict,
            CloudType& cloud,
            const List<label>& associatedPatches
        );


    //- Destructor
    virtual ~Current();


    virtual autoPtr<BoundaryEvent<CloudType>> clone() const
    {
        return autoPtr<BoundaryEvent<CloudType>>
        (
            new Current<CloudType>(*this)
        );
    }

    // Member Functions

        //- Check deletion
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

        //- Check insertion
        virtual void ejectionEvent(typename CloudType::parcelType& p, label patchId);

        //- Print statistics
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Current.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
