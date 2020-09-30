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

\*---------------------------------------------------------------------------*/

#include "picParcel.H"
#include "PICCloud.H"
#include "FreeStream.H"
#include "NoInflow.H"
#include "ZARMInOutflow.H"
#include "SimpleInjection.H"
#include "IdealCurrentSource.H"
#include "FloatingPotential.H"
#include "ParticleEmitter.H"
#include "circuitRLC.H"
#include "ThermionicEmission.H"
#include "NoEmission.H"
#include "FowlerNordheim.H"
#include "IdealVoltageSource.H"
#include "SputterEmission.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef PICCloud<picParcel> CloudType;

    //Define ParticleEmitter
    template class ParticleEmitter<CloudType>;

    //Emission models
    makeEmissionModel(CloudType)
    makeEmissionModelType(ThermionicEmission,CloudType)
    makeEmissionModelType(FowlerNordheim,CloudType)
    makeEmissionModelType(SputterEmission,CloudType)
    makeEmissionModelType(NoEmission,CloudType)

    //Open/Freestream Boundaries
    makeBoundaryModel(CloudType)
    makeBoundaryModelType(NoInflow, CloudType)
    makeBoundaryModelType(FreeStream, CloudType)
    makeBoundaryModelType(ZARMInOutflow, CloudType)

    //Simple particle injection
    makeBoundaryModelType(SimpleInjection, CloudType)

    //Circuit Boundaries
    makeBoundaryModelType(circuitRLC, CloudType) //General Circuit
    makeBoundaryModelType(FloatingPotential, CloudType) //Open Circuit
    makeBoundaryModelType(IdealCurrentSource, CloudType)//Current-Driven (Open Circuit)
    makeBoundaryModelType(IdealVoltageSource, CloudType)//Short Circuit
}


// ************************************************************************* //
