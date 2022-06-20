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
    Foam::CoulombCollisionModel

Description
    Templated coulomb collision class

SourceFiles
    CoulombCollisionModel.C
    CoulombCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CoulombCollisionModel_H
#define CoulombCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//Pairing-Submodel
template<class CloudType>
class PairingAlgorithm;

template<class CloudType>
class IonizationModel;

template<class CloudType>
class WeightCorrectionModel;

/*---------------------------------------------------------------------------*\
                     Class CoulombCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CoulombCollisionModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Pointer to PairingAlgorithm
        autoPtr<PairingAlgorithm<CloudType>>
            pairingAlgorithm_;

        //- Pointer to IonizationModel
        autoPtr<IonizationModel<CloudType>>
            ionizationModel_;

        //- Allow collisions between particles of one species?
        bool allowIntraCollision_;

        //- Calculate the debye length or are the coulomb logs given?
        bool calculateDebyeLength_;

        //- Pointer to the WeightCorrection model
        autoPtr<WeightCorrectionModel<CloudType>>
            weightCorrection_;

protected:
    // Protected data

        //- Storage of the coulomb logarithm for all species combinations
        List<List<scalar>> coulombLog_; //there are n(n+1)/2 pairs but I think using two list n*n entries results in faster lookup: list[s1][s2]

        //- Calculate the average Coulomb log per time step using this variable
        scalar average_coulombLog_;

        //- Account for drift in the temperature calculation used for the calculation of the Debye length (default: true)
        Switch debyeLength_accountForDrift_;

        //- Only consider these species in the calculation of the Debye length
        List<label> debyeLength_Species_;

        //- Calculate the average impact parameter per time step using this variable
        scalar average_b_;

        //- Consider b_min in the calculation of the impact parameter (default: true)
        Switch calculate_bmin_;
public:

    //- Runtime type information
    TypeName("CoulombCollisionModel")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CoulombCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    )


    // Constructors

        //- Construct null from owner
        CoulombCollisionModel(CloudType& owner);

        //- Construct from components
        CoulombCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~CoulombCollisionModel();


    //- Selector
    static autoPtr<CoulombCollisionModel<CloudType>> New
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

        //- Return allowIntraCollision_
        inline const bool& allowIntraCollision() const;

        //- Return calculateDebyeLength_
        inline const bool& calculateDebyeLength() const;

        //- Return reference to PairingAlgorithm
        inline const PairingAlgorithm<CloudType>&
            pairingAlgorithm() const;

        //- Return reference to IonizationModel
        inline const IonizationModel<CloudType>&
            ionizationModel() const;

        //- Return non-const reference to IonizationModel
        inline IonizationModel<CloudType>&
            ionizationModel();

        //- Return reference to WeightCorrectionModel
        inline WeightCorrectionModel<CloudType>&
            weightCorrection();

        //- Return non-const reference to WeightCorrectionModel
        inline const WeightCorrectionModel<CloudType>&
            weightCorrection() const;

        //- Return non-const reference to the average Coulomb log. Used to update the value in child classes.
        inline scalar& average_coulombLog();

        //- Return const reference to average Coulomb log
        inline const scalar& average_coulombLog() const;

        //- Return non-const reference to average impact parameter b
        inline scalar& average_impactParameter();

        //- Return const reference to average impact parameter b
        inline const scalar& average_impactParameter() const;

    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const = 0;

        //- Calculate and return the debye length in celli
        virtual scalar debyeLength(label celli);

        //- Apply collision
        virtual void collide
        (
            typename CloudType::parcelType* pP,
            typename CloudType::parcelType* pQ,
            scalar debyeLength,
            scalar nP,
            scalar nQ,
            scalar nPQ
        ) = 0;

        //- Pair and collide all charged particle
        void handleCollisions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "CoulombCollisionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCoulombCollisionModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CoulombCollisionModel<CloudType>,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        CoulombCollisionModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeCoulombCollisionModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    CoulombCollisionModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CoulombCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
