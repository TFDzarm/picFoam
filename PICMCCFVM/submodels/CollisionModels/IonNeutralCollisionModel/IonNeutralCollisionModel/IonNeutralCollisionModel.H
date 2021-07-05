/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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
    Foam::IonNeutralCollisionModel

Description
    Templated PIC particle ion neutral collision class

SourceFiles
    IonNeutralCollisionModel.C
    IonNeutralCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef IonNeutralCollisionModel_H
#define IonNeutralCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "CrossSectionModelList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IonNeutralCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WeightCorrectionModel;

template<class CloudType>
class IonNeutralCollisionModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- PointerList to CrossSectionModels
        CrossSectionList<CloudType, crossSectionType::IonElasticCS> elasticCSModels_;
        CrossSectionList<CloudType, crossSectionType::IonChargeExCS> chargeExCSModels_;

        //- Maximum total cross section *cR for background gas collisions
        volScalarField backgroundSigmaTcRMax_;

        //-  Remainder for parcel selection (background)
        volScalarField::Internal backgroundCollisionRemainder_;

        //- Maximum total cross section *cR
        volScalarField sigmaTcRMax_;

        //-  Remainder for parcel selection FIXME: only works for one neutral species
        volScalarField::Internal collisionRemainder_;

        //- Pointer to WeightCorrection model, used for unequal parcel weight
        autoPtr<WeightCorrectionModel<CloudType>>
            weightCorrection_;

        //- Should we update the velocity of the neutral parcel? (def: true)
        bool updateNeutralVelocity_;
public:

    //- Runtime type information
    TypeName("IonNeutralCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        IonNeutralCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        IonNeutralCollisionModel(CloudType& owner);

        //- Construct from components
        IonNeutralCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~IonNeutralCollisionModel();


    //- Selector
    static autoPtr<IonNeutralCollisionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
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

        //- Return PointerList of elastic cross section models
        inline const CrossSectionList<CloudType,crossSectionType::IonElasticCS>&
            elasticCrossSections() const;

        //- Return PointerList of excitation cross section models
        inline const CrossSectionList<CloudType,crossSectionType::IonChargeExCS>&
            chargeExCrossSections() const;

        //- Return wheter to update the neutrals velocity during collision or not
        inline const bool& updateNeutralVelocity() const;

    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const = 0;      

        //- Elastic collision
        void elasticCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ);

        //- Elastic background collision
        void elasticCollision(typename CloudType::parcelType* pP, vector& Uq);

        //- Charge exchange collision
        void chargeExchangeCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ);

        //- Charge exchange background collision
        void chargeExchangeCollision(typename CloudType::parcelType* pP, vector& Uq);

        //- Update the velocity through the choosen method
        virtual void updateVelocity
        (
            typename CloudType::parcelType& pP,
            typename CloudType::parcelType& pQ
        ) = 0;

        //- Update the velocity through the choosen method (for background collisions)
        virtual void updateVelocity
        (
            typename CloudType::parcelType& pP,
            vector& Uq,
            label idQ
        ) = 0;


        //- Select particles and collide
        void handleCollisions();

        //- Initialize
        void initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities);

private:
        //- Collision process
        void performCollisions();

        //- Background collision process
        void performBackgroundCollisions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IonNeutralCollisionModelI.H"

#define makeIonNeutralCollisionModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        IonNeutralCollisionModel<CloudType>,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        IonNeutralCollisionModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeIonNeutralCollisionModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    IonNeutralCollisionModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IonNeutralCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
