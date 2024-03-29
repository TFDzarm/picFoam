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
#include "CloudCompositionInfo.H"
#include "MomentumInfo.H"
#include "KineticEnergyInfo.H"
#include "TemperatureInfo.H"
#include "FieldEnergyInfo.H"
#include "IonizationInfo.H"
#include "PrintParcelInfo.H"
#include "DumpInfo.H"
#include "VelocityInfo.H"
#include "UnidirectionalVelocityInfo.H"
#include "ChargeInfo.H"
#include "CellOccupancyInfo.H"
#include "IsotropizationInfo.H"
#include "DebyeLengthInfo.H"
#include "PlasmaFrequencyInfo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef PICCloud<picParcel> CloudType;

    //Define the selection table
    makeDiagnosticInfo(CloudType);

    //Add the Diagnostic models
    makeDiagnosticInfoType(CloudCompositionInfo, CloudType);
    makeDiagnosticInfoType(MomentumInfo, CloudType);
    makeDiagnosticInfoType(KineticEnergyInfo, CloudType);
    makeDiagnosticInfoType(TemperatureInfo, CloudType);
    makeDiagnosticInfoType(FieldEnergyInfo, CloudType);
    makeDiagnosticInfoType(IonizationInfo, CloudType);
    makeDiagnosticInfoType(PrintParcelInfo, CloudType);
    makeDiagnosticInfoType(DumpInfo, CloudType);
    makeDiagnosticInfoType(VelocityInfo, CloudType);
    makeDiagnosticInfoType(UnidirectionalVelocityInfo, CloudType);
    makeDiagnosticInfoType(ChargeInfo, CloudType);
    makeDiagnosticInfoType(CellOccupancyInfo, CloudType);
    makeDiagnosticInfoType(IsotropizationInfo, CloudType);
    makeDiagnosticInfoType(DebyeLengthInfo, CloudType);
    makeDiagnosticInfoType(PlasmaFrequencyInfo, CloudType);
}

// ************************************************************************* //
