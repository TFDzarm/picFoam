/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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

#include "VelocityInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VelocityInfo<CloudType>::VelocityInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    vDrift_(cloud.typeIdList().size(),Zero),
    nParticleTypes_(cloud.typeIdList().size(),0.0)
{}

template<class CloudType>
Foam::VelocityInfo<CloudType>::VelocityInfo(const VelocityInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    vDrift_(im.vDrift_),
    nParticleTypes_(im.nParticleTypes_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VelocityInfo<CloudType>::~VelocityInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VelocityInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    vDrift_[p.typeId()] += p.U()*p.nParticle();
    nParticleTypes_[p.typeId()] += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::VelocityInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    Pstream::listCombineGather(vDrift_, plusEqOp<vector>());
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());

    Info << "    Drift Velocities" << nl;
    forAll(cloud.typeIdList(),id)
    {
        if(nParticleTypes_[id] > 0)
        {
            vector vDrift = vDrift_[id]/nParticleTypes_[id];
            Info << "    [" << cloud.typeIdList()[id] << "]                             = " << vDrift << " == |" << mag(vDrift) << "| m/s" << nl;
        }
    }

   //Reset list
   forAll(cloud.typeIdList(),i)
   {
       vDrift_[i] = Zero;
       nParticleTypes_[i] = 0.0;
   }
}

// ************************************************************************* //
