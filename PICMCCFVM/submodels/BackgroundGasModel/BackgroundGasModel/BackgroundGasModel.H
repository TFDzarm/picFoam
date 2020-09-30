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
    Foam::BackgroundGasModel

Description
    Templated background gas model class

SourceFiles
    BackgroundGasModel.C
    BackgroundGasModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef BackgroundGasModel_H
#define BackgroundGasModel_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BackgroundGasModel Declaration
\*---------------------------------------------------------------------------*/

enum class CollisionEvent { Elastic, Exciation, Ionization, ChargeExchange };

template<class CloudType>
class BackgroundGasModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("BackgroundGasModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BackgroundGasModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        BackgroundGasModel(CloudType& owner);

        //- Construct from components
        BackgroundGasModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~BackgroundGasModel();


    //- Selector
    static autoPtr<BackgroundGasModel<CloudType>> New
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


    // Member Functions

        //- Is model active?
        virtual bool active() const = 0;

        //- Update the background gas, called before collision()
        virtual void update() = 0;

        //- Called when a parcel collied with the background
        virtual void collisionUpdate(CollisionEvent type, label celli, const vector& preU, const vector& postU) = 0;

        //- The background species
        virtual label species() const = 0;

        //- The number of equivalent particles used in ionizations
        virtual scalar nParticle() const = 0;

        //- The background gas fields
        virtual tmp<scalarField> numberDensity() const = 0;
        virtual tmp<scalarField> temperature() const = 0;

        //-  Velocity of a background gas parcel used during collision
        virtual vector sampleVelocity(label celli) = 0;

        //- Initalize
        virtual void initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBackgroundGasModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(BackgroundGasModel<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        BackgroundGasModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeBackgroundGasModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    BackgroundGasModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BackgroundGasModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
