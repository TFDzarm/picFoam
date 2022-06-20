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

#include "IsotropizationInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropizationInfo<CloudType>::IsotropizationInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    vSqr_(cloud.typeIdList().size(),Zero),
    vDrift_(cloud.typeIdList().size(),Zero),
    nParticleTypes_(cloud.typeIdList().size(),0.0),
    accountForDrift_(true)//default
{
    accountForDrift_.readIfPresent("accountForDrift",this->coeffDict());
    if(accountForDrift_)
        Info << "       |= Do not incorporate dirft velocities" << endl;
    else
        Info << "       |= Incorporate the drift velocities in calculations" << endl;
}

template<class CloudType>
Foam::IsotropizationInfo<CloudType>::IsotropizationInfo(const IsotropizationInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    vSqr_(im.vSqr_),
    vDrift_(im.vDrift_),
    nParticleTypes_(im.nParticleTypes_),
    accountForDrift_(im.accountForDrift_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropizationInfo<CloudType>::~IsotropizationInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IsotropizationInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const vector& U(p.U());
    //Calculate all required properties
    vSqr_[p.typeId()] += cmptMultiply(U,U)*p.nParticle();
    vDrift_[p.typeId()] += U*p.nParticle();
    nParticleTypes_[p.typeId()] += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::IsotropizationInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());
    
    //Parallel COM the properties
    Pstream::listCombineGather(vSqr_, plusEqOp<vector>());
    Pstream::listCombineGather(vDrift_, plusEqOp<vector>());
    Pstream::listCombineGather(nParticleTypes_,plusEqOp<scalar>());

    //Info on every species ... only includes drift if chosen   
    forAll(cloud.typeIdList(),id)
    {
        if(nParticleTypes_[id] > 0)
        {
            vector vDrift = vDrift_[id]/nParticleTypes_[id];
            scalar vDrift2x = 0;//default
            scalar vDrift2yz = 0;//default

            if(accountForDrift_)//Do we add the drift velocity?
            {
                vDrift2x = vDrift.x()*vDrift.x();
                vDrift2yz = vDrift.y()*vDrift.y()+vDrift.z()*vDrift.z();
            }

            scalar Tx = cloud.constProps(id).mass()/(constant::physicoChemical::k.value())*(vSqr_[id].x()/nParticleTypes_[id] - vDrift2x);
            scalar Tyz = cloud.constProps(id).mass()/(2.0*constant::physicoChemical::k.value())*((vSqr_[id].y()+vSqr_[id].z())/nParticleTypes_[id] - vDrift2yz);
            scalar f = Tyz/Tx;
            Info << "    [" << cloud.typeIdList()[id] << "]" << nl
                 << "    T_x                             = " << Tx << " K == " << Tx*constant::physicoChemical::k.value()/constant::electromagnetic::e.value() << " eV" << nl
                 << "    T_yz                            = " << Tyz << " K == " << Tyz*constant::physicoChemical::k.value()/constant::electromagnetic::e.value() << " eV" << nl
                 << "    fraction                        = " << f << nl;

        }
    }

   //Reset list
   forAll(cloud.typeIdList(),i)
   {
       vSqr_[i] = Zero;
       vDrift_[i] = Zero;
       nParticleTypes_[i] = 0.0;
   }
}

// ************************************************************************* //
