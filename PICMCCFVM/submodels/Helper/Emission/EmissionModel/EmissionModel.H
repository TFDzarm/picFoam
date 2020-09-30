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
    Foam::EmissionModel

Description
    Template EmissionModel class

SourceFiles
    EmissionModel.C
    EmissionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef EmissionModel_H
#define EmissionModel_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "ParticleEmitter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EmissionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EmissionModel: public ParticleEmitter<CloudType>
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

public:

    //- Runtime type information
    TypeName("EmissionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        EmissionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        EmissionModel(CloudType& cloud);

        //- Construct from components
        EmissionModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

    //- Construct and return a clone
    virtual autoPtr<EmissionModel<CloudType>> clone() const = 0;

    //- Destructor
    virtual ~EmissionModel();


    //- Selector
    static autoPtr<EmissionModel<CloudType>> New
    (
        const word& modelType,
        const dictionary& dict,
        CloudType& cloud
    );


    //- Return the dictionary
    const dictionary& dict() const;

    //- Return the coefficients dictionary
    const dictionary& coeffDict() const;


    // Member Functions

        //- Emission of particles from boundary
        virtual void emission() = 0;

        //- Collisional emission e.g. sputter
        virtual void collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeEmissionModel(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(EmissionModel<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        EmissionModel<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeEmissionModelType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    EmissionModel<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EmissionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
