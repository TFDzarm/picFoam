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
    Foam::RajuElasticCS

Description
    Elastic cross section model for argon.

    G.G. Raju
    Electron-atom Collision Cross Sections in Argon: An Analysis and Comments
    https://doi.org/10.1109/TDEI.2004.1324355

\*---------------------------------------------------------------------------*/

#ifndef RajuElasticCS_H
#define RajuElasticCS_H

#include "CrossSectionModel.H"
#include "ListListOps.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class RajuElasticCS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, crossSectionType Type>
class RajuElasticCS
:
    public CrossSectionModel<CloudType,Type>
{
//Private data

    //- Storage for interpolated cross sections
    List<scalar> elasticCS_;

    //- Energies provided by Raju
    List<scalar> raju_Energies = {0.08,0.1,0.2,0.5,1,5,10,12,14,15,16,18,20,25,30,40,50,70,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000};

    //- Total cross section data...unused
    //List<scalar> raju_sigmaQt = {1.5,1.19,0.5,0.374,1.07,8.81,20.06,21.31,22.797,23.474,22.963,20.433,19.738,16.45,16.362,
    //                             12.414,10.73,9.327,8.237,8.054,6.748,5.83,5.207,4.622,4.277,3.968,3.702,3.479,3.306,3.151,2.988,2.843,2.701,2.603,2.512,2.429,2.268};

    //- Elastic cross sections provided by Raju
    List<scalar> raju_Qel = {1.5,1.19,0.5,0.374,1.07,8.81,20.06,21.25,22.56,23.22,22.3,19.52,18.6,16.27, \
                             13.94,9.51,7.74,6.15,5.0,4.86,3.79,3.2,2.833,2.466,2.29,2.115,2.0,1.885,1.806,1.728,1.649,1.571,1.492,1.448,1.405,1.361,1.274};
public:

    //- Runtime type information
    TypeName("Raju");


    // Constructors

        //- Construct from dictionary
        RajuElasticCS
        (
            const dictionary& dict,
            CloudType& cloud,
            const label& associatedTypeId
        );


    //- Destructor
    virtual ~RajuElasticCS();

    virtual autoPtr<CrossSectionModel<CloudType,Type>> clone() const
    {
            return autoPtr<CrossSectionModel<CloudType,Type>>
            (
                new RajuElasticCS<CloudType,Type>(*this)
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
    #include "RajuElastic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
