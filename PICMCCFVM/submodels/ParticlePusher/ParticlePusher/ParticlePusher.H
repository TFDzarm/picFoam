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
    Foam::ParticlePusher

Description
    Templated particle pusher model class

SourceFiles
    ParticlePusher.C
    ParticlePusherNew.C

\*---------------------------------------------------------------------------*/

#ifndef ParticlePusher_H
#define ParticlePusher_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ParticlePusher Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticlePusher
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
    TypeName("ParticlePusher");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticlePusher,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        ParticlePusher(CloudType& owner);

        //- Construct from components
        ParticlePusher
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~ParticlePusher();


    //- Selector
    static autoPtr<ParticlePusher<CloudType>> New
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

        //- Update the parcels velocity
        virtual void updateVelocity(typename CloudType::parcelType& p, const scalar trackTime) = 0;

        //- Return the velocity change
        virtual vector correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParticlePusher(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(ParticlePusher<CloudType>, 0);   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        ParticlePusher<CloudType>,                                       \
        dictionary                                                             \
    );


#define makeParticlePusherType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    ParticlePusher<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticlePusher.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
