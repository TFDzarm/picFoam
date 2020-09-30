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
    Foam::EmissionModelList

Description
    List of emission models

SourceFiles
    EmissionModelList.C

\*---------------------------------------------------------------------------*/

#ifndef EmissionModelList_H
#define EmissionModelList_H

#include "PtrList.H"
#include "EmissionModel.H"
#include "ParticleEmitter.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class InjectionModelList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EmissionModelList
:
    public PtrList<EmissionModel<CloudType>>
{
public:

    // Constructors

        //- Construct null from owner
        EmissionModelList(CloudType& owner);

        //- Construct from dictionary and cloud owner
        EmissionModelList(const dictionary& dict, CloudType& owner);


        //- Construct copy
        EmissionModelList(const EmissionModelList<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<EmissionModelList<CloudType>> clone() const
        {
            return autoPtr<EmissionModelList<CloudType>>
            (
                new EmissionModelList<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~EmissionModelList();

    //Member functions

        //- Setup models
        void setupModels(const dictionary& dict, CloudType& owner);

        //- Call initilizeParticleEmitter() of the ParticleEmitter
        void initilizeAll(label patchId, typename ParticleEmitter<CloudType>::VelocityModel model);

        //- Call emission() on all models
        void emission();

        //- Call collisionalEmission() on all models
        void collisionalEmission(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EmissionModelList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


