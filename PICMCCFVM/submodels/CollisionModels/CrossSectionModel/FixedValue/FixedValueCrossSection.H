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
    Foam::FixedValueCrossSection

Description
    Simple fixed value cross section used for testing.

\*---------------------------------------------------------------------------*/

#ifndef FixedValueCrossSection_H
#define FixedValueCrossSection_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class FixedValueCrossSection Declaration
\*---------------------------------------------------------------------------*/

//Base Template
template<class CloudType, crossSectionType Type>
class FixedValueCrossSection
:
    public CrossSectionModel<CloudType,Type>
{
public:

    //- Runtime type information
    TypeName("FixedValue");


    // Constructors

        //- Construct from dictionary
        FixedValueCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~FixedValueCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new FixedValueCrossSection<CloudType,Type>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};

//Template specialization Elastic
template <class CloudType>
class FixedValueCrossSection<CloudType,crossSectionType::ElectronElasticCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>
{
//Private data

    //- the fixed value
    scalar value_;
public:
    //- Runtime type information
    TypeName("FixedValue");


    // Constructors

        //- Construct from dictionary
        FixedValueCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~FixedValueCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>>
            (
                new FixedValueCrossSection<CloudType,crossSectionType::ElectronElasticCS>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;


};


//Template specialization Excitation
template <class CloudType>
class FixedValueCrossSection<CloudType,crossSectionType::ElectronExciationCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>
{
//Private data

    //- the fixed value and threshold
    scalar value_;
    scalar threshold_;
public:
    //- Runtime type information
    TypeName("FixedValue");


    // Constructors

        //- Construct from dictionary
        FixedValueCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~FixedValueCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>>
            (
                new FixedValueCrossSection<CloudType,crossSectionType::ElectronExciationCS>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;


};


//Template specialization Ionization
template <class CloudType>
class FixedValueCrossSection<CloudType,crossSectionType::ElectronIonizationCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>
{
//Private data

    //- the fixed value and threshold
    scalar value_;
    scalar threshold_;
public:
    //- Runtime type information
    TypeName("FixedValue");


    // Constructors

        //- Construct from dictionary
        FixedValueCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~FixedValueCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>>
            (
                new FixedValueCrossSection<CloudType,crossSectionType::ElectronIonizationCS>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedValueCrossSection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
