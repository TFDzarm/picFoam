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
    Foam::BoundaryModel


Description
    Templated BoundaryModel class

SourceFiles
    BoundaryModel.C
    BoundaryModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef BoundaryModel_H
#define BoundaryModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BoundaryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoundaryModel
{
protected:
    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        // Reference to the owner cloud class
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
    TypeName("BoundaryModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BoundaryModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const List<label>& associatedPatches
        ),
        (dict, owner, associatedPatches)
    );


    // Constructors

        //- Construct null from owner
        BoundaryModel(CloudType& owner, const List<label>& associatedPatches);

        BoundaryModel
        (
            const dictionary& dict,
            CloudType& owner,
            const List<label>& associatedPatches
        );

        //- Construct from dictionary
        BoundaryModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type,
            const List<label>& associatedPatches
        );

        //- Construct and return a clone
        virtual autoPtr<BoundaryModel<CloudType>> clone() const = 0;

    //- Destructor
    virtual ~BoundaryModel();


    //- Selector
    static autoPtr<BoundaryModel<CloudType>> New
    (
        const word& modelType,
        const dictionary& dict,
        CloudType& owner,
        const List<label>& associatedPatches
    );


    inline bool interactWithPatch(label patchId);

    // Access

        //- Return const access the owner cloud object
        inline const CloudType& owner() const;

        //- Return non-const access the owner cloud object for manipulation
        inline CloudType& owner();

        //- Return the owner cloud dictionary
        inline const dictionary& dict() const;

        //- Return the coefficients dictionary
        inline const dictionary& coeffDict() const;



    //Member functions

        //- Introduce particles
        virtual void injection() = 0;

        //- Called before and after Maxwell solver
        virtual void preUpdate_Boundary() { }
        virtual void postUpdate_Boundary() { }

        //- Patch boundary condition (e.g. reflect, delete, ...), return true if interaction occurred
        virtual bool particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
        {
            return false;
        }

        //- Called if particle was ejected on the patch (also called by BoundaryEvent)
        virtual void particleEjection(typename CloudType::parcelType& p, label patchId) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "BoundaryModelI.H"

#define makeBoundaryModel(CloudType)                                     \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(BoundaryModel<CloudType>, 0);    \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        BoundaryModel<CloudType>,                                        \
        dictionary                                                             \
    );


#define makeBoundaryModelType(SS, CloudType)                             \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    BoundaryModel<CloudType>::                                           \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BoundaryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
