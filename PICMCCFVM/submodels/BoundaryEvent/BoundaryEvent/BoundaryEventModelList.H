/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::BoundaryEventModelList

Description
    List of boundary event models

SourceFiles
    BoundaryEventModelList.C

\*---------------------------------------------------------------------------*/

#ifndef BoundaryEventModelList_H
#define BoundaryEventModelList_H

#include "PtrList.H"
#include "BoundaryEvent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BoundaryEventModelList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoundaryEventModelList
:
    public PtrList<BoundaryEvent<CloudType>>
{
public:

    // Constructors

        //- Construct null from owner
        BoundaryEventModelList(CloudType& owner);

        //- Construct from dictionary and cloud owner
        BoundaryEventModelList(const dictionary& dict, CloudType& owner);

        //- Construct copy
        BoundaryEventModelList(const BoundaryEventModelList<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<BoundaryEventModelList<CloudType>> clone() const
        {
            return autoPtr<BoundaryEventModelList<CloudType>>
            (
                new BoundaryEventModelList<CloudType>(*this)
            );
        }

        //- Setup all models stored in the list
        void setupModels(const dictionary& dict, CloudType& owner);

    //- Destructor
    virtual ~BoundaryEventModelList();


    // Member Functions

        //- Call collisionEvent on all associated models
        void onCollision
        (
            typename CloudType::parcelType& p,
            typename CloudType::parcelType::trackingData& td
        );

        //- Call ejectionEvent on all associated models
        void onEjection
        (
            typename CloudType::parcelType& p,
            label patchId
        );

        //- Call postMove on all models
        void postMove();

        //- Call info on all models
        void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BoundaryEventModelList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


