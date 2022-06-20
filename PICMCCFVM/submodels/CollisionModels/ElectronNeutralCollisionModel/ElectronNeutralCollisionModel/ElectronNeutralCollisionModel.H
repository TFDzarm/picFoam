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
    Foam::ElectronNeutralCollisionModel

Description
    Templated PIC particle electron neutral collision class

SourceFiles
    ElectronNeutralCollisionModel.C
    ElectronNeutralCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ElectronNeutralCollisionModel_H
#define ElectronNeutralCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "CrossSectionModelList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BinaryElectronCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WeightCorrectionModel;

template<class CloudType>
class ElectronNeutralCollisionModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- PointerList to CrossSectionModels
        CrossSectionList<CloudType, crossSectionType::ElectronElasticCS> elasticCSModels_;
        CrossSectionList<CloudType, crossSectionType::ElectronExciationCS> excitationCSModels_;
        CrossSectionList<CloudType, crossSectionType::ElectronIonizationCS> ionizationCSModels_;

        //- Should we delte neutral after ionization? (def: true)
        bool deleteNeutral_;

        //- Should we update the velocity of the neutral parcel? (def: true)
        bool updateNeutralVelocity_;

        scalar ionizationNEquivalentParticles_;

        //- Maximum total cross section *cR for background gas collisions
        volScalarField backgroundSigmaTcRMax_;

        //-  Remainder for parcel selection (background)
        volScalarField::Internal backgroundCollisionRemainder_;


        //- Maximum total cross section *cR per neutral species
        PtrList<volScalarField> sigmaTcRMax_;

        //-  Remainder for parcel selection per neutral species
        PtrList<volScalarField::Internal> collisionRemainder_;

        //- Pointer to WeightCorrection model, used for unequal parcel weight
        autoPtr<WeightCorrectionModel<CloudType>>
            weightCorrection_;

        //- Running average of the ionization rate
        scalar ionizationAverage_;

        //- Time when the averaging process began
        scalar averageStart_;

        //- Fields to keep track to where these event occur
        volScalarField elasticCollisions_;
        volScalarField excitationCollisions_;
        volScalarField ionizationCollisions_;

        //- Do we create the new ion or electron during ionization? (def: true)
        bool createElectron_;
        bool createIon_;
public:

    //- Runtime type information
    TypeName("ElectronNeutralCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ElectronNeutralCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        ElectronNeutralCollisionModel(CloudType& owner);

        //- Construct from components
        ElectronNeutralCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~ElectronNeutralCollisionModel();


    //- Selector
    static autoPtr<ElectronNeutralCollisionModel<CloudType>> New
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

        //- Return whether to delete neutrals after ionization or not
        inline const bool& deleteNeutrals() const;

        //- Return whether to update the neutrals velocity during collision or not
        inline const bool& updateNeutralVelocity() const;

        //- Return whether to create the new ion/electron or not
        inline const bool& createIon() const;
        inline const bool& createElectron() const;

        //- Return PointerList of elastic cross section models
        inline const CrossSectionList<CloudType,crossSectionType::ElectronElasticCS>&
            elasticCrossSections() const;

        //- Return PointerList of excitation cross section models
        inline const CrossSectionList<CloudType,crossSectionType::ElectronExciationCS>&
            excitationCrossSections() const;

        //- Return PointerList of ionization cross section models
        inline const CrossSectionList<CloudType,crossSectionType::ElectronIonizationCS>&
            ionizationCrossSections() const;

    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const = 0;

        //- Handle excitation collision event
        void excitationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN);
        void excitationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, vector& Un);

        //- Handle ionization collision event
        bool ionizationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN);
        bool ionizationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, vector& Un);

        //- Handle elastic collision event
        void elasticCollision(scalar eV, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN);
        void elasticCollision(scalar eV, typename CloudType::parcelType* pE, vector& Un);

        //- Update the velocities during collisions with the selected model
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, typename CloudType::parcelType& pN) = 0;
        virtual void updateVelocity(scalar eV, typename CloudType::parcelType& pE, vector& Un, label idN) = 0;

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

#include "ElectronNeutralCollisionModelI.H"

#define makeElectronNeutralCollisionModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        ElectronNeutralCollisionModel<CloudType>,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        ElectronNeutralCollisionModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeElectronNeutralCollisionModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    ElectronNeutralCollisionModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ElectronNeutralCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
