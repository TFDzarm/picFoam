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

#include "MomentumInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MomentumInfo<CloudType>::MomentumInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    momentumofTypes_(cloud.typeIdList().size(),Zero)
{}

template<class CloudType>
Foam::MomentumInfo<CloudType>::MomentumInfo(const MomentumInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    momentumofTypes_(im.momentumofTypes_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MomentumInfo<CloudType>::~MomentumInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MomentumInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const CloudType& cloud(this->owner());
    const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(p.typeId());

    momentumofTypes_[p.typeId()] += cP.mass()*p.U()*p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::MomentumInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    vector linearMomentum = cloud.linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());
    Pstream::listCombineGather(momentumofTypes_, plusEqOp<vector>());

    Info << "   |Total linear momentum|          = "
         << mag(linearMomentum) << nl;

   forAll(cloud.typeIdList(),id)
   {
       vector momentum = momentumofTypes_[id];
       Info<< "    [" << cloud.typeIdList()[id] << "]" << "                             = "
           << mag(momentum) << " Ns" << nl;
   }

   //Reset list
   forAll(momentumofTypes_,i)
       momentumofTypes_[i] = vector::zero;

}

// ************************************************************************* //
