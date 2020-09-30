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
    Foam::BoundaryModelList

Description
    List of boundary models

SourceFiles
    BoundaryModelList.C

\*---------------------------------------------------------------------------*/

#ifndef BoundaryModelList_H
#define BoundaryModelList_H

#include "PtrList.H"
#include "BoundaryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BoundaryModelList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoundaryModelList
:
    public PtrList<BoundaryModel<CloudType>>
{
public:

    // Constructors

        //- Construct null from owner
        BoundaryModelList(CloudType& owner);

        //- Construct from dictionary and cloud owner
        BoundaryModelList(const dictionary& dict, CloudType& owner);

        //- Construct copy
        BoundaryModelList(const BoundaryModelList<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<BoundaryModelList<CloudType>> clone() const
        {
            return autoPtr<BoundaryModelList<CloudType>>
            (
                new BoundaryModelList<CloudType>(*this)
            );
        }

        //- Setup all models stored in the list
        void setupModels(const dictionary& dict, CloudType& owner);

    //- Destructor
    virtual ~BoundaryModelList();



    // Member Functions

       //- Call injection on all models
       void injection();

       //- Call preUpdate_Boundary on all models
       void preUpdate_Boundary();

       //- Call postUpdate_Boundary on all models
       void postUpdate_Boundary();

       //- Call particleBC on all models
       bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);

       //- Call particleEjection on all models
       void onEjection(typename CloudType::parcelType& p, label patchId);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BoundaryModelList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


