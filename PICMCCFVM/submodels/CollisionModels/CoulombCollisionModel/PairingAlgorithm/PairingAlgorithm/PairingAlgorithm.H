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
    Foam::PairingAlgorithm

Description
    Templated coulomb collision pairing class

SourceFiles
    PairingAlgorithm.C
    PairingAlgorithmNew.C

\*---------------------------------------------------------------------------*/

#ifndef PairingAlgorithm_H
#define PairingAlgorithm_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class CloudType>
class CoulombCollisionModel;

/*---------------------------------------------------------------------------*\
                     Class PairingAlgorithm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PairingAlgorithm
{
    // Private data

        //- The CoulombCollisionModel dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& cloud_;

        //- Reference to the CoulombCollisionModel
        Foam::CoulombCollisionModel<CloudType>& coulombCollisionModel_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Keep track of the global cell count, we require it for the Debye length average
        label nCells_;
public:

    //- Runtime type information
    TypeName("PairingAlgorithm")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PairingAlgorithm,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud,
            Foam::CoulombCollisionModel<CloudType>& collisionModel
        ),
        (dict, cloud, collisionModel)
    )


    // Constructors

        //- Construct null from owner
        PairingAlgorithm(CloudType& cloud, Foam::CoulombCollisionModel<CloudType>& collisionModel);

        //- Construct from components
        PairingAlgorithm
        (
            const dictionary& dict,
            CloudType& cloud,
            Foam::CoulombCollisionModel<CloudType>& collisionModel,
            const word& type
        );


    //- Destructor
    virtual ~PairingAlgorithm();


    //- Selector
    static autoPtr<PairingAlgorithm<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner,
        Foam::CoulombCollisionModel<CloudType>& collisionModel
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

        //- Return the owner CoulombCollision object
        inline const Foam::CoulombCollisionModel<CloudType>& coulombCollisionModel() const;

        //- Return non-const access to the owner CoulombCollision object
        inline Foam::CoulombCollisionModel<CloudType>& coulombCollisionModel();

        //- Return the global cell count
        inline const label& nCellsGlobal() const;

    // Member Functions

        //- Select particle pairs and collide
        virtual void pairANDcollide() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "PairingAlgorithmI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePairingAlgorithm(CloudType)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        PairingAlgorithm<CloudType>,                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        PairingAlgorithm<CloudType>,                                       \
        dictionary                                                             \
    );


#define makePairingAlgorithmType(SS, CloudType)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    PairingAlgorithm<CloudType>::                                          \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PairingAlgorithm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
