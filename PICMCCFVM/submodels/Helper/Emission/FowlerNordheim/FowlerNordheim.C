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

#include "FowlerNordheim.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FowlerNordheim<CloudType>::FowlerNordheim
(
    const dictionary& dict,
    CloudType& cloud
)
:
    EmissionModel<CloudType>(dict,cloud,typeName),
    A_(1.541434e-6),//A eV V^-2
    B_(6.830889776e9),//eV^-3/2 V m^-1
    t_(1.0),//1.1
    fnW_(readScalar(this->coeffDict().lookup("fnW"))),
    nParticle_(this->coeffDict().lookupOrDefault("nParticle",-1.0)),
    pRemainder_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FowlerNordheim<CloudType>::~FowlerNordheim()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FowlerNordheim<CloudType>::emission()
{
    CloudType& cloud(this->cloud());
    const fvMesh& mesh(cloud.mesh());

    const scalar& electronTypeId(cloud.electronTypeId());

    scalar electronNequiv;
    if(nParticle_ > 0.0)
        electronNequiv = nParticle_;
    else
        electronNequiv = cloud.constProps(electronTypeId).nParticle();

    const scalar& electronMagCharge = mag(cloud.constProps(electronTypeId).charge());

    scalar dt = mesh.time().deltaTValue();

    scalar j = 0.0;
    const volVectorField::Boundary& EBf(cloud.electricField().boundaryFieldRef());
    forAll(EBf[this->patchId()],lF)//calculate the current density on every face of the patch
    {
        vector Ei = EBf[this->patchId()][lF];
        if((Ei&this->patchNormals()[lF]) > 0.0)//patchNormal points out of the domain
        {
            //Fowler Nordheim eq. with Wang, Loew approxcimation
            //important at > 0.5 GV/m with fnW = 2
            scalar magE = mag(Ei);
            scalar y=3.7947e-5*sqrt(magE)/fnW_;
            scalar v=0.956-1.062*y*y;

            j += A_*magE*magE/fnW_/t_/t_*exp(-B_*v*sqrt(fnW_)*sqrt(fnW_)*sqrt(fnW_)/magE);
        }
    }
    //Parallel COM sum up j_
    reduce(j,sumOp<scalar>());

    //Calculate number of particles to be injected
    scalar nEmission = j*this->patchArea()*dt/(electronNequiv*electronMagCharge) + pRemainder_;
    label nEmit = label(nEmission);
    pRemainder_ = nEmission-scalar(nEmit);

    //Print some info
    Info << "    FowlerNordheim emission " << nEmit << " parcels (rmdr: " << pRemainder_ << ")" << endl;

    //Emit particles call into emitter class
    for(label i=0; i<nEmit; i++)
    {
        this->emitParticle(fnW_*electronMagCharge/constant::physicoChemical::k.value(),electronTypeId,electronNequiv);
    }
}


// ************************************************************************* //
