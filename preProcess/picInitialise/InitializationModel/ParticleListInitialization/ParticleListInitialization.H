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
    Foam::ParticleListInitialization

Description
    Initalize particle from a list

\*---------------------------------------------------------------------------*/

#ifndef ParticleListInitialization_H
#define ParticleListInitialization_H

#include "InitializationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class ParticleListInitialization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleListInitialization
:
    public InitializationModel<CloudType>
{
//Private data

    //- Is temperature given in electronvolt?
    bool useElectronVolt_;
public:

    //- Runtime type information
    TypeName("ParticleList");


    // Constructors

        //- Construct from dictionary
        ParticleListInitialization
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~ParticleListInitialization();


    // Member Functions

        //- Initalize...
        virtual void initialiseParticles();

        //- Set by model do we read before or calculate after initialization
        virtual bool calculateNumberDensities()
        {
            return true;
        }

        //- Set by model do we read before or calculate after initialization
        virtual bool calculateTemperatures()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleListInitialization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

