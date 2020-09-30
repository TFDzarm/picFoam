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
    Foam::BoundaryEvent

Description
    Templated BoundaryEvent model class

SourceFiles
    BoundaryEvent.C
    BoundaryEventNew.C

\*---------------------------------------------------------------------------*/

#ifndef BoundaryEvent_H
#define BoundaryEvent_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BoundaryEvent Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoundaryEvent
{
protected:
    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;
private:
    // Private data

        //- List of patchIds that use this model
        List<label> associatedPatches_;

        //- List of all patches (convenience function, true if patch uses the model)
        List<bool> interactionTable_;

        //- Build interactionTable_ from associatedPatches_
        void initInteractionTable();
public:

    //- Runtime type information
    TypeName("BoundaryEvent")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BoundaryEvent,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const List<label>& associatedPatches
        ),
        (dict, owner, associatedPatches)
    )


    // Constructors

        //- Construct null from owner
        BoundaryEvent(CloudType& owner);


        BoundaryEvent
        (
            const dictionary& dict,
            CloudType& owner,
            const List<label>& associatedPatches
        );

        //- Construct from components
        BoundaryEvent
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type,
            const List<label>& associatedPatches
        );

        //- Construct and return a clone
        virtual autoPtr<BoundaryEvent<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~BoundaryEvent();


    //- Selector
    static autoPtr<BoundaryEvent<CloudType>> New
    (
        const word& modelType,
        const dictionary& dict,
        CloudType& owner,
        const List<label>& associatedPatches

    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return non-const access to the owner cloud object
        CloudType& owner();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Return associatedPatches_
        const List<label>& associatedPatches() const;

        inline bool interactWithPatch(label patchId);

    // Member Functions

        //- Called when parcel collides with an associated patch
        virtual void collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td) = 0;

        //- Called if a parcel is ejected from an associated patch
        virtual void ejectionEvent(typename CloudType::parcelType& p, label patchId) {}

        //- Called after all parcels have been moved
        virtual void postMove() {}

        //- Called after timestep finished
        virtual void info() {} 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoundaryEventI.H"

#define makeBoundaryEvent(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(BoundaryEvent<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        BoundaryEvent<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeBoundaryEventType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    BoundaryEvent<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BoundaryEvent.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
