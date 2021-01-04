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

#include "ThermionicEmission.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermionicEmission<CloudType>::ThermionicEmission
(
    const dictionary& dict,
    CloudType& cloud
)
:
    EmissionModel<CloudType>(dict,cloud,typeName),
    lambdaR_(readScalar(this->coeffDict().lookup("lambdaR"))),
    A0_(1201732.101),//4*pi*m*k^2*q/h^3
    Tbc_(readScalar(this->coeffDict().lookup("T"))),
    fnW_(readScalar(this->coeffDict().lookup("fnW"))),
    j_(0.0),
    nParticle_(this->coeffDict().lookupOrDefault("nParticle",-1.0)),
    pRemainder_(0.0)
{
    //Calculate the current density for the given temperature
    const scalar& electronTypeId(cloud.electronTypeId());
    const scalar& electronMagCharge = mag(cloud.constProps(electronTypeId).charge());
    j_ = lambdaR_*A0_*Tbc_*Tbc_*exp(-fnW_*electronMagCharge/(Tbc_*constant::physicoChemical::k.value()));

    Info << "|->    J = " << j_ << " A/m^2" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermionicEmission<CloudType>::~ThermionicEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermionicEmission<CloudType>::emission()
{
    CloudType& cloud(this->cloud());
    const fvMesh& mesh(cloud.mesh());

    const scalar& electronTypeId(cloud.electronTypeId());
    scalar electronNequiv;
    if(nParticle_ > 0.0)//If we specified the weight assign it
        electronNequiv = nParticle_;
    else
        electronNequiv = cloud.constProps(electronTypeId).nParticle();

    const scalar& electronMagCharge = mag(cloud.constProps(electronTypeId).charge());

    scalar dt = mesh.time().deltaTValue();

    //Calculate number of particles to be emitted
    scalar nEmission = j_*this->patchArea()*dt/(electronNequiv*electronMagCharge) + pRemainder_;
    label nEmit = label(nEmission);
    pRemainder_ = nEmission-scalar(nEmit);

    Info << "    Thermionic emission " << nEmit << " parcels (rmdr: " << pRemainder_ << ")" << endl;
    for(label i=0; i<nEmit; i++)
    {
        this->emitParticle(Tbc_,electronTypeId,electronNequiv);//Creates electron and updates boundary models
    }
}


// ************************************************************************* //
