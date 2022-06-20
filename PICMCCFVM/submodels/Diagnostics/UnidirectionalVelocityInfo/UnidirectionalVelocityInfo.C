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

#include "UnidirectionalVelocityInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UnidirectionalVelocityInfo<CloudType>::UnidirectionalVelocityInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    directionMask_(this->coeffDict().lookup("direction")),
    velocitySqr_(cloud.typeIdList().size(),Zero),
    nParticleTypes_(cloud.typeIdList().size(),0.0)
{
    //Create mask used for the dot product with the velocity vector thus only selected directions are includes in the magnitude
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if(directionMask_[cmpt] != 0.0)
            directionMask_[cmpt] = 1.0;
    }
    Info << "       |= Calculate velocity magnitude in direction " << directionMask_ << endl;
}

template<class CloudType>
Foam::UnidirectionalVelocityInfo<CloudType>::UnidirectionalVelocityInfo(const UnidirectionalVelocityInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    directionMask_(im.directionMask_),
    velocitySqr_(im.velocitySqr_),
    nParticleTypes_(im.nParticleTypes_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UnidirectionalVelocityInfo<CloudType>::~UnidirectionalVelocityInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::UnidirectionalVelocityInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    nParticleTypes_[p.typeId()] += p.nParticle();

    //Keep the vector, thus make a component multiply here
    velocitySqr_[p.typeId()] += p.nParticle()*cmptMultiply(p.U(),p.U());
}

//- Print info
template<class CloudType>
void Foam::UnidirectionalVelocityInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    //Parallel communication
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());
    Pstream::listCombineGather(velocitySqr_,plusEqOp<vector>());

    Info << "    <v^2> in direction " << directionMask_ << nl;
    forAll(cloud.typeIdList(),id)
    {
        if(nParticleTypes_[id] > 0)
        {
            //Do the dot product to consider only user selected directions
            scalar v = (velocitySqr_[id]&directionMask_)/nParticleTypes_[id];
            Info << "    [" << cloud.typeIdList()[id] << "]                             = " << v << " m/s" << nl;
        }
    }

   //Reset list
   forAll(cloud.typeIdList(),i)
   {
       velocitySqr_[i] = Zero;
       nParticleTypes_[i] = 0.0;
   }
}

// ************************************************************************* //
