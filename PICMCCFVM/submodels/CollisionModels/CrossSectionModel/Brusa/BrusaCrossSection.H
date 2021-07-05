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
    Foam::BrusaCrossSectionBase

Description
    Elastic cross section for Neon, Argon, Krypton, Xenon.

    R.S. Brusa et al.
    Analytical partitioning of total cross sections for electron scattering on noble gases
    https://doi.org/10.1007/s004600050092

\*---------------------------------------------------------------------------*/

#ifndef BrusaCrossSection_H
#define BrusaCrossSection_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class BrusaCrossSectionBase Declaration
\*---------------------------------------------------------------------------*/

//Base class
class BrusaCrossSectionBase
{
protected:
//Protected data

    //- Species that can be selected
    enum species
    {
        undefined = -1,
        Ne,Ar,Kr,Xe
    };

    //- Selected species
    species s_;
public:
    // Constructors
        //- Construct from dictionary
        BrusaCrossSectionBase
        (
            const dictionary& dict
        );

        virtual ~BrusaCrossSectionBase();
};

//Base Template
template<class CloudType, crossSectionType Type>
class BrusaCrossSection
:
    public BrusaCrossSectionBase, public CrossSectionModel<CloudType,Type>
{
public:

    //- Runtime type information
    TypeName("Brusa");


    // Constructors

        //- Construct from dictionary
        BrusaCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BrusaCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new BrusaCrossSection<CloudType,Type>(*this)
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
class BrusaCrossSection<CloudType,crossSectionType::ElectronElasticCS>
:
    public BrusaCrossSectionBase, public CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>
{
//Private data

    //- Data provided by Brusa
    scalarList R1_ = {7.076,0.632,4.115,1.210};
    scalarList S1_ = {0.0381,0.593,0.0017,0.203};
    scalarList R2_ = {3.64,3.65,0.388,0.121};
    scalarList S2_ = {1.11,0,1.35,13.6};

    scalarList minEel_ = {10.0,10.0,10.0,80.0};
    scalarList maxEel_ = {10000.0,10000.0,10000.0,10000.0};

public:
    //- Runtime type information
    TypeName("Brusa");


    // Constructors

        //- Construct from dictionary
        BrusaCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BrusaCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronElasticCS>>
            (
                new BrusaCrossSection<CloudType,crossSectionType::ElectronElasticCS>(*this)
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
class BrusaCrossSection<CloudType,crossSectionType::ElectronExciationCS>
:
    public BrusaCrossSectionBase, public CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>
{
//Private data

    //- Data provided by Brusa
    scalarList F_ = {77.65,85.97,25.19,22.0,18.27};
    scalarList G_ = {0,31.7e-3,23.6e-3,23.3e-3,0};
    scalarList Eex_ = {19.8,16.619,11.5,9.915,8.315};
    scalarList minEex_ = {20.0,30.0,20.0,20.0,80.0};
    scalarList maxEex_ = {500.0,4000.0,3000.0,4000.0,10000.0};

public:
    //- Runtime type information
    TypeName("Brusa");


    // Constructors

        //- Construct from dictionary
        BrusaCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BrusaCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronExciationCS>>
            (
                new BrusaCrossSection<CloudType,crossSectionType::ElectronExciationCS>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};


//TemplateSpecialization Ionization
template <class CloudType>
class BrusaCrossSection<CloudType,crossSectionType::ElectronIonizationCS>
:
    public BrusaCrossSectionBase, public CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>
{
//Private data

    //- Data provided by Brusa
    scalarList L_ = {2000.0,7.92,78.76,33.76,1000.0};
    scalarList M_ = {70.9,7.04,18.62,15.93,53.79};
    scalarList N_ = {0.0,0.0,25.66,14.45,109.6};
    scalarList P_ = {1.49e-3,21.6e-3,8.42e-3,12.6e-3,3.58e-3};
    scalarList Eion_ = {24.587,21.584,15.759,13.999,12.130};
    scalarList minEion_ = {24.587,21.584,15.759,13.999,12.130};
    scalarList maxEion_ = {5000.0,5000.0,5000.0,5000.0,5000.0};

public:
    //- Runtime type information
    TypeName("Brusa");


    // Constructors

        //- Construct from dictionary
        BrusaCrossSection
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BrusaCrossSection();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::ElectronIonizationCS>>
            (
                new BrusaCrossSection<CloudType,crossSectionType::ElectronIonizationCS>(*this)
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
    #include "BrusaCrossSection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
