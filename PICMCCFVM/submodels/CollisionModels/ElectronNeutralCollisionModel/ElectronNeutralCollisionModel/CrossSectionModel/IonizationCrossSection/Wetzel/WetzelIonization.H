/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::HardSphere

Description
    Hard Sphere BinaryCollision Model

\*---------------------------------------------------------------------------*/

#ifndef WetzelIonizationCS_H
#define WetzelIonizationCS_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class HardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class WetzelIonizationCs
:
    public CrossSectionModel<CloudType,Type>
{
    List<scalar> ionizationCS_;

    List<scalar> wetzel_Energies = {15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0,120.0,125.0,130.0,135.0,140.0,145.0,150.0,155.0,160.0,165.0,170.0,175.0,180.0,185.0,190.0,195.0,200.0};

    List<scalar> wetzel_Qi1 = {0.0,0.03,0.15,0.29,0.44,0.58,0.71,0.85,0.96,1.1,1.25,1.34,1.44,1.57,1.68,1.74,1.86,1.92,2.02,2.07,2.2,2.24,2.31,2.39,2.42,2.47,2.56,2.62,2.62,2.63,2.63,2.65,2.66,2.68,2.7,2.68,2.67,2.68,2.67,2.64,2.61,2.59,2.57,2.56,2.54,2.54,2.53,2.5,2.48,2.47,2.44,2.45,2.42,2.4,2.39,2.32};

public:

    //- Runtime type information
    TypeName("Wetzel");


    // Constructors

        //- Construct from dictionary
        WetzelIonizationCs
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~WetzelIonizationCs();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new WetzelIonizationCs<CloudType,Type>(*this)
            );
    }


    // Member Functions
    virtual scalar crossSection(scalar eVEnergy) const;

    virtual scalar threshold() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WetzelIonization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
