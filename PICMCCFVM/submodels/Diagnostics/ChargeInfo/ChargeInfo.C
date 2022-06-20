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

#include "ChargeInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeInfo<CloudType>::ChargeInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName)
{}

template<class CloudType>
Foam::ChargeInfo<CloudType>::ChargeInfo(const ChargeInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeInfo<CloudType>::~ChargeInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ChargeInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{}

//- Print info
template<class CloudType>
void Foam::ChargeInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    tmp<scalarField> totalChargeDensity(cloud.rhoCharge()+cloud.spaceChargeDensity());
    scalar totalCharge = sum(totalChargeDensity/cloud.mesh().cellVolumes());

    reduce(totalCharge, sumOp<scalar>());

    Info << "   Total charge                     = " << totalCharge << " As" << nl;
}

// ************************************************************************* //
