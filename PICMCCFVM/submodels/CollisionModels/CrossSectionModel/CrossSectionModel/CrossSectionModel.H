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
    Foam::CrossSectionModel

Description
    Templated wall cross section model class. For electron neutral collisions.

SourceFiles
    CrossSectionModel.C
    CrossSectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CrossSectionModel_H
#define CrossSectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CrossSectionModel Declaration
\*---------------------------------------------------------------------------*/

enum class crossSectionType
{
    ElectronElasticCS,
    ElectronExciationCS,
    ElectronIonizationCS,
    IonElasticCS,
    IonChargeExCS
};

template<class CloudType, crossSectionType Type>
class CrossSectionModel
{
    // Private data

    //- The cloud dictionary
    const dictionary& dict_;

    // reference to the owner cloud class
    CloudType& owner_;

    //- The coefficients dictionary
    const dictionary coeffDict_;

    //- Associated typeId
    const label associatedTypeId_;

public:
    //- Runtime type information
    TypeName("CrossSectionModel")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CrossSectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const label& associatedTypeId
        ),
        (dict, owner, associatedTypeId)
    )

    // Constructors

    //- Construct null from owner
    CrossSectionModel(CloudType& owner, const label& associatedTypeId);


    CrossSectionModel
    (
        const dictionary& dict,
        CloudType& owner,
        const label& associatedTypeId
    );

    //- Construct from components
    CrossSectionModel
    (
        const dictionary& dict,
        CloudType& owner,
        const word& type,
        const label& associatedTypeId
    );

    virtual autoPtr<CrossSectionModel<CloudType, Type>> clone() const = 0;

    //- Destructor
    virtual ~CrossSectionModel();

    //- Selector
    static autoPtr<CrossSectionModel<CloudType, Type>> New
    (
        const word& modelName,
        const dictionary& dict,
        CloudType& owner,
        const label associatedTypeId

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

    //- Return the associated typeid
    label associatedTypeId() const;

    // Member Functions

    //- Calculate the cross section
    virtual scalar crossSection(scalar eVEnergy) const = 0;

    //- Return the threshold
    virtual scalar threshold() const = 0;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeCSModel(CrossSectionType)                                    \
                                   \
    defineNamedTemplateTypeNameAndDebug(CrossSectionType, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        CrossSectionType,                                       \
        dictionary                                                             \
    );


#define makeCSType(SS, CloudType, Type)                            \
    typedef SS<CloudType,Type> cs##SS##Template;                                                                           \
    defineNamedTemplateTypeNameAndDebug(cs##SS##Template, 0);                     \
                                                                               \
    CrossSectionModel<CloudType,Type>::                                          \
        adddictionaryConstructorToTable<cs##SS##Template>                         \
            add##SS##ConstructorToTable_;

//We need this when SS is the same classname for multiple events, the string ##SS## wouldnt be unique in usage of makeCSType
#define makeCSTypeMulti(SS, CrossSectionType, CloudType, Type)                            \
    typedef SS<CloudType,Type> cs##SS##CrossSectionType##Template;                                                                           \
    defineNamedTemplateTypeNameAndDebug(cs##SS##CrossSectionType##Template, 0);                     \
                                                                               \
    CrossSectionModel<CloudType,Type>::                                          \
        adddictionaryConstructorToTable<cs##SS##CrossSectionType##Template>                         \
            add##SS##CrossSectionType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CrossSectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
