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
    Foam::ChargeDistribution

Description
    Templated ChargeDistribution class

SourceFiles
    ChargeDistribution.C
    ChargeDistributionNew.C

\*---------------------------------------------------------------------------*/

#ifndef ChargeDistribution_H
#define ChargeDistribution_H

#include "fvCFD.H"
#include "barycentric.H"
#include "tetIndices.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ChargeDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ChargeDistribution
{
// Private data

    //- The cloud dictionary
    const dictionary& dict_;

    // reference to the owner cloud class
    CloudType& cloud_;

    //- The coefficients dictionary
    const dictionary coeffDict_;

protected:
// Protected data

    //- The charge density field
    volScalarField field_;

public:

    //- Runtime type information
    TypeName("ChargeDistribution");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ChargeDistribution,
        dictionary,
        (
            const word& fieldName,
            const dictionary& dict,
            CloudType& cloud
        ),
        (fieldName, dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        ChargeDistribution(const word& fieldName, CloudType& owner);

        //- Construct from components
        ChargeDistribution
        (
            const word& fieldName,
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~ChargeDistribution();


    //- Selector
    static autoPtr<ChargeDistribution<CloudType>> New
    (
        const word& fieldName,
        const dictionary& dict,
        CloudType& cloud
    );


    // Access

        //- Return the owner cloud object
        const CloudType& cloud() const;

        //- Return non-const access to the owner cloud object
        CloudType& cloud();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Return the charge density field
        inline volScalarField& field();
        inline const volScalarField& field() const;

    // Member Functions

        //- Add particle charge to the field
        virtual void add(const typename CloudType::parcelType& parcel) = 0;

        //- Reset the charge density field
        virtual void reset();

        //- Update the distribution
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChargeDistributionI.H"

#define makeChargeDistribution(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(ChargeDistribution<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        ChargeDistribution<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeChargeDistributionType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    ChargeDistribution<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChargeDistribution.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
