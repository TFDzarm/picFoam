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
    Foam::StraubIonizationCS

Description
    Ionization cross section model for argon.

    H.C. Straub et al.
    Absolute partial and total cross sections for electron-impact ionization of argon from threshold to 1000 eV
    https://doi.org/10.1103/PhysRevA.52.1115

\*---------------------------------------------------------------------------*/

#ifndef StraubIonizationCS_H
#define StraubIonizationCS_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class HardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class StraubIonizationCS
:
    public CrossSectionModel<CloudType,Type>
{
//Private data

    //- Storage for interpolated cross sections
    List<scalar> ionizationCS_;

    //- Energies provided by Straub
    List<scalar> straub_Energies = {17.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,110.0,120.0,140.0,160.0,180.0,200.0,225.0,250.0,275.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.0};

    //- Ionization cross sections provided by Straub
    List<scalar> straub_Qi1 = {0.017,0.46,1.24,1.84,2.26,2.55,2.66,2.70,2.69,2.67,2.67,2.67,2.66,2.69,2.70,2.69,2.67,2.64,2.61,2.55,2.45,2.35,2.27,2.18,2.10,1.99,1.87,1.79,1.63,1.53,1.39,1.31,1.23,1.16,1.09,1.03,0.976,0.932,0.901,0.865,0.824,0.795};

public:

    //- Runtime type information
    TypeName("Straub");


    // Constructors

        //- Construct from dictionary
        StraubIonizationCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~StraubIonizationCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new StraubIonizationCS<CloudType,Type>(*this)
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
    #include "StraubIonization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
