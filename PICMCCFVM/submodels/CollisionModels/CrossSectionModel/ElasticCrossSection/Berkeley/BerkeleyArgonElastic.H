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
    Foam::BerkeleyArgonElasticCS

Description
    Elastic cross section model for Argon.
    Copied from XPDP1: https://ptsg.egr.msu.edu

\*---------------------------------------------------------------------------*/

#ifndef BerkeleyArgonElastic_H
#define BerkeleyArgonElastic_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class HardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class BerkeleyArgonElasticCS
:
    public CrossSectionModel<CloudType,Type>
{
//Private Data

    //- This model uses pre-calculated values for interpolation
    List<scalar> elasticLowEnergyCS_;
    List<scalar> elasticCS_;

public:

    //- Runtime type information
    TypeName("BerkeleyArgon");


    // Constructors

        //- Construct from dictionary
        BerkeleyArgonElasticCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BerkeleyArgonElasticCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new BerkeleyArgonElasticCS<CloudType,Type>(*this)
            );
    }


    // Member Functions

        //- Calculate the cross section
        virtual scalar crossSection(scalar eVEnergy) const;

        //- Return the threshold
        virtual scalar threshold() const;

};


//Ion-Neutral Collision
template<class CloudType>
class BerkeleyArgonElasticCS<CloudType,crossSectionType::IonElasticCS>
:
    public CrossSectionModel<CloudType,crossSectionType::IonElasticCS>
{
//Private Data
public:

    //- Runtime type information
    TypeName("BerkeleyArgon");


    // Constructors

        //- Construct from dictionary
        BerkeleyArgonElasticCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~BerkeleyArgonElasticCS();

    virtual autoPtr<CrossSectionModel<CloudType,crossSectionType::IonElasticCS>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,crossSectionType::IonElasticCS>>
            (
                new BerkeleyArgonElasticCS<CloudType,crossSectionType::IonElasticCS>(*this)
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
    #include "BerkeleyArgonElastic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
