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
    Foam::RajuExcitationCS

Description
    Excitation cross section model for argon.

    G.G. Raju
    Electron-atom Collision Cross Sections in Argon: An Analysis and Comments
    https://doi.org/10.1109/TDEI.2004.1324355

\*---------------------------------------------------------------------------*/

#ifndef RajuExcitationCS_H
#define RajuExcitationCS_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class RajuExcitationCS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class RajuExcitationCS
:
    public CrossSectionModel<CloudType,Type>
{
//Private data

    //- Storage for interpolated cross sections
    List<scalar> excitationCS_;

    //- Energies provided by Raju
    List<scalar> raju_Energies = {0.08,0.1,0.2,0.5,1,5,10,12,14,15,16,18,20,25,30,40,50,70,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000};

    //- Excitation cross sections provided by Raju
    List<scalar> raju_Qex = {0.045,0.237,0.254,0.262,0.565,0.643,0.676,0.619,0.511,0.457,0.406,0.378,0.344,0.275,0.237,0.201,0.177, \
                                0.175,0.173,0.154,0.134,0.128,0.121,0.116,0.111,0.101,0.091,0.087,0.083,0.079};
public:

    //- Runtime type information
    TypeName("Raju");


    // Constructors

        //- Construct from dictionary
        RajuExcitationCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~RajuExcitationCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new RajuExcitationCS<CloudType,Type>(*this)
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
    #include "RajuExcitation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
