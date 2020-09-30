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
    Foam::BinaryCollisionModel

Description
    Templated PIC particle collision class. For neutral and ion species.

SourceFiles
    BinaryCollisionModel.C
    BinaryCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef BinaryCollisionModel_H
#define BinaryCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CloudType>
class TotalCrossSectionModel;

template<class CloudType>
class WeightCorrectionModel;
/*---------------------------------------------------------------------------*\
                     Class BinaryCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BinaryCollisionModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- Separate group for ions
        volScalarField sigmaTIoncRMax_;       

        //- A field holding the remainder from the previous collision selections
        volScalarField::Internal collisionSelectionRemainder_;

        //- A field holding the remainder from the previous collision selections for light parcels
        volScalarField::Internal collisionSelectionIonRemainder_;

        //- Pointer to the TotalCorssSection model
        autoPtr<TotalCrossSectionModel<CloudType>>
            totalCrossSection_;

        //- Propbability of an ion to collide with the background gas model
        volScalarField::Internal backgroundCollisionPropIon_;

        //- Remainder of ion - background collision candidates
        volScalarField::Internal backgroundCSRIon_;

        //- Propbability an of neutral to collide with the background gas model
        volScalarField::Internal backgroundCollisionProp_;

        //- Remainder of neutral - background collision candidates
        volScalarField::Internal backgroundCSR_;

        //- Pointer to the WeightCorrection model used for unequal particle weights
        autoPtr<WeightCorrectionModel<CloudType>>
            weightCorrection_;

        //- Switch whether to use the simple charge exchange model by Nanbu: https://doi.org/10.1109/27.887765
        bool chargeExchangeCollision_;
public:

    //- Runtime type information
    TypeName("BinaryCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BinaryCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        BinaryCollisionModel(CloudType& owner);

        //- Construct from components
        BinaryCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~BinaryCollisionModel();


    //- Selector
    static autoPtr<BinaryCollisionModel<CloudType>> New
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

        //- Return the sigmaTcRMax field.  non-const access to allow
        // updating.
        inline volScalarField& sigmaTcRMax();

        //- Return the sigmaTcRMax field.  non-const access to allow
        // updating.
        inline volScalarField& sigmaTIoncRMax();

        //- Return the collision selection remainder field.  non-const
        // access to allow updating.
        inline scalarField& collisionSelectionRemainder();

        //- Return the ion collision selection remainder field.  non-const
        // access to allow updating.
        inline scalarField& collisionSelectionIonRemainder();

    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const = 0;

        //- Elastic collision
        void elasticCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ);

        //- Elastic background collision
        void elasticCollision(typename CloudType::parcelType* pP, vector& Uq, label idQ);

        //- Charge exchange collision
        void chargeExchangeCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ);

        //- Charge exchange background collision
        void chargeExchangeCollision(typename CloudType::parcelType* pP, vector& Uq, label idQ);

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

        //- Start collision process if active
        void handleCollisions();

        //- Init sigmaTcRMax
        void initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities);
private:
        //- Choose particles and collide
        void performCollisions();

        //- Choose particles and collide with background
        void performBackgroundCollisions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "BinaryCollisionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBinaryCollisionModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        BinaryCollisionModel<CloudType>,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        BinaryCollisionModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeBinaryCollisionModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    BinaryCollisionModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BinaryCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
