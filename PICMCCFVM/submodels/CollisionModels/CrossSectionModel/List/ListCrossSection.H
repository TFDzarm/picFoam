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
    Foam::ListCrossSection

Description
    Simple fixed value cross section used for testing.

\*---------------------------------------------------------------------------*/

#ifndef ListCrossSection_H
#define ListCrossSection_H

#include "CrossSectionModel.H"
#include "Field.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class ListCrossSection Declaration
\*---------------------------------------------------------------------------*/

//Base Template
template<class CloudType, crossSectionType Type>
class ListCrossSection
:
    public CrossSectionModel<CloudType,Type>
{
public:

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Construct from dictionary
        ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new ListCrossSection<CloudType,Type>(*this)
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
class ListCrossSection<CloudType,crossSectionType::ElectronElasticCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>
{
//Private data

    //- the values
    Field<scalar> energies_;
    Field<scalar> values_;

public:
    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Construct from dictionary
        ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>>
            (
                new ListCrossSection<CloudType,crossSectionType::ElectronElasticCS>(*this)
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
class ListCrossSection<CloudType,crossSectionType::ElectronExciationCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>
{
//Private data

    //- the fixed value and threshold
    Field<scalar> energies_;
    Field<scalar> values_;
    scalar threshold_;
public:
    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Construct from dictionary
        ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>>
            (
                new ListCrossSection<CloudType,crossSectionType::ElectronExciationCS>(*this)
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
class ListCrossSection<CloudType,crossSectionType::ElectronIonizationCS>
:
    public CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>
{
//Private data

    //- the fixed value and threshold
    Field<scalar> energies_;
    Field<scalar> values_;
    scalar threshold_;
public:
    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Construct from dictionary
        ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>>
            (
                new ListCrossSection<CloudType,crossSectionType::ElectronIonizationCS>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};



//Template specialization IonElastic
template <class CloudType>
class ListCrossSection<CloudType,crossSectionType::IonElasticCS>
    :
    public CrossSectionModel<CloudType,crossSectionType::IonElasticCS>
{
    //Private data

    //- the fixed value
    Field<scalar> energies_;
    Field<scalar> values_;
public:
    //- Runtime type information
    TypeName("List");


    // Constructors

    //- Construct from dictionary
    ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
            );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::IonElasticCS>> clone() const
    {
        return autoPtr<CrossSectionModel<CloudType,crossSectionType::IonElasticCS>>
            (
                new ListCrossSection<CloudType,crossSectionType::IonElasticCS>(*this)
                );
    }


    // Member Functions

    //- Calculate the cross section
    virtual scalar crossSection(scalar eVEnergy) const;

    //- Return the threshold
    virtual scalar threshold() const;

};


//Template specialization IonChargeExchange
template <class CloudType>
class ListCrossSection<CloudType,crossSectionType::IonChargeExCS>
    :
    public CrossSectionModel<CloudType,crossSectionType::IonChargeExCS>
{
    //Private data

    //- the fixed value and threshold
    Field<scalar> energies_;
    Field<scalar> values_;
    scalar threshold_;
public:
    //- Runtime type information
    TypeName("List");


    // Constructors

    //- Construct from dictionary
    ListCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
            );


    //- Destructor
    virtual ~ListCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::IonChargeExCS>> clone() const
    {
        return autoPtr<CrossSectionModel<CloudType,crossSectionType::IonChargeExCS>>
            (
                new ListCrossSection<CloudType,crossSectionType::IonChargeExCS>(*this)
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
    #include "ListCrossSection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
