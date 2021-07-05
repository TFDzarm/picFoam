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
    Foam::RajuIonizationCS

Description
    Ionization cross section model for argon.

    G.G. Raju
    Electron-atom Collision Cross Sections in Argon: An Analysis and Comments
    https://doi.org/10.1109/TDEI.2004.1324355

\*---------------------------------------------------------------------------*/

#ifndef RajuIonizationCS_H
#define RajuIonizationCS_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class HardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class RajuIonizationCS
:
    public CrossSectionModel<CloudType,Type>
{
//Private data

    //- Storage for interpolated cross sections
    List<scalar> ionizationCS_;

    //- Energies provided by Raju
    List<scalar> raju_Energies = {0.08,0.1,0.2,0.5,1,5,10,12,14,15,16,18,20,25,30,40,50,70,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000};

    //- Ionization cross sections provided by Raju
    List<scalar> raju_Qi = {0.02,0.294,0.627,1.302,1.803,2.393,2.533,2.771,2.859,2.85,2.683,2.393,2.173,1.979,1.812,1.68,1.548, \
                                1.46,1.372,1.302,1.223,1.161,1.108,1.064,1.02,0.985,0.915};

public:

    //- Runtime type information
    TypeName("Raju");


    // Constructors

        //- Construct from dictionary
        RajuIonizationCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~RajuIonizationCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new RajuIonizationCS<CloudType,Type>(*this)
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
    #include "RajuIonization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
