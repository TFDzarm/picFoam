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
    Foam::InitializationModel

Description
    Template initialization model class.

\*---------------------------------------------------------------------------*/

#ifndef InitializationModel_H
#define InitializationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

namespace Foam
{

template<class CloudType>
class InitializationModel
{
    //- The initialise dictionary
    const dictionary& dict_;

    //- Reference to the owner cloud class
    CloudType& cloud_;

    //- The coefficients dictionary
    const dictionary coeffDict_;

    //- Index of the cell zone we wish to initalize in
    label cellZoneIndex_;
protected:

    //- Number densities (set or calculated)
    Field<scalar> numberDensities_;

    //- Temperatures (set or calculated)
    Field<scalar> temperatures_;

public:
    InitializationModel
    (
        const dictionary& dict,
        CloudType& owner
    );

    InitializationModel
    (
        const dictionary& dict,
        CloudType& owner,
        const word& type
    );

    //- Destructor
    virtual ~InitializationModel();

    TypeName("InitializationModel")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InitializationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    )

    //- Selector
    static autoPtr<InitializationModel<CloudType>> New
    (
        const word& modelType,
        const dictionary& dict,
        CloudType& cloud
    );

    // Access

        //- Return const access the owner cloud object
        inline const CloudType& cloud() const;

        //- Return non-const access the owner cloud object for manipulation
        inline CloudType& cloud();

        //- Return the owner cloud dictionary
        inline const dictionary& dict() const;

        //- Return the coefficients dictionary
        inline const dictionary& coeffDict() const;


        //- Set by model do we read before or calculate after initialization
        virtual bool calculateNumberDensities()
        {
            return false;
        }

        //- Set by model do we read before or calculate after initialization
        virtual bool calculateTemperatures()
        {
            return false;
        }

        //- Return number density
        inline Field<scalar>& numberDensities();

        //- Return const number densities
        inline const Field<scalar>& numberDensities() const;

        //- Return temperature
        inline Field<scalar>& temperatures();

        //- Return const temperature
        inline const Field<scalar>& temperatures() const;

        //- Return cell zone id
        inline const label& cellZone() const;

        //- Function to read number density and temperature if needed
        void readDefaultSettings();


        //- Initialise...
        virtual void initialiseParticles() = 0;
};

} // End namespace Foam

#define makeInitializationModel(CloudType)                                     \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(InitializationModel<CloudType>, 0);    \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        InitializationModel<CloudType>,                                        \
        dictionary                                                             \
    );


#define makeInitializationModelType(SS, CloudType)                             \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    InitializationModel<CloudType>::                                           \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


#ifdef NoRepository
    #include "InitializationModel.C"
#endif

#endif
