/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
#include "EquipartitionInitialization.H"
#include "QuasineutralPlasmaInitialization.H"
#include "ParticleListInitialization.H"
#include "NoParticleInitialization.H"
#include "IsotropizationInitialization.H"
#include "CenterInitialization.H"
#include "BeamInitialization.H"
#include "CrossSectionInitialization.H"
#include "CountInitialization.H"
#include "UniformInitialization.H"
#include "MaxwellianInitialization.H"

namespace Foam
{
    typedef PICCloud<picParcel> CloudType;

    makeInitializationModel(CloudType);

    // Add instances of the initialization model to the table
    makeInitializationModelType(IsotropizationInitialization, CloudType);
    makeInitializationModelType(EquipartitionInitialization, CloudType);
    makeInitializationModelType(QuasineutralPlasmaInitialization, CloudType);
    makeInitializationModelType(ParticleListInitialization, CloudType);
    makeInitializationModelType(NoParticleInitialization, CloudType);
    makeInitializationModelType(CenterInitialization, CloudType);
    makeInitializationModelType(BeamInitialization, CloudType);
    makeInitializationModelType(CrossSectionInitialization, CloudType);
    makeInitializationModelType(CountInitialization, CloudType);
    makeInitializationModelType(UniformInitialization, CloudType);
    makeInitializationModelType(MaxwellianInitialization, CloudType);
}
