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

#include "PlasmaFrequencyInfo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlasmaFrequencyInfo<CloudType>::PlasmaFrequencyInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    nParticle_(0.0),
    totalVolume_(sum(cloud.mesh().cellVolumes())),
    warnTimeStep_(false)
{
    const polyMesh& mesh(cloud.mesh());
    totalVolume_ = sum(mesh.cellVolumes());
    reduce(totalVolume_,sumOp<scalar>());

    if(cloud.electronTypeId() < 0)
        FatalErrorInFunction << "No electron species defined." << abort(FatalError);

    warnTimeStep_.readIfPresent("warnTimeStep",this->coeffDict());
    if(warnTimeStep_)
        Info << "       |= Warn when the time step is too larger compared to the plasma frequency" << endl;
}

template<class CloudType>
Foam::PlasmaFrequencyInfo<CloudType>::PlasmaFrequencyInfo(const PlasmaFrequencyInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    nParticle_(im.nParticle_),
    totalVolume_(im.totalVolume_),
    warnTimeStep_(im.warnTimeStep_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PlasmaFrequencyInfo<CloudType>::~PlasmaFrequencyInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PlasmaFrequencyInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{
    const CloudType& cloud(this->owner());

    //Sum the number of electron particles
    if(p.typeId() == cloud.electronTypeId())
        nParticle_ += p.nParticle();
}

//- Print info
template<class CloudType>
void Foam::PlasmaFrequencyInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());
    const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(cloud.electronTypeId());

    //Parallel COM
    reduce(nParticle_,sumOp<scalar>());

    scalar n = nParticle_/totalVolume_;
    scalar plasmaFrequency = sqrt(n*cP.charge()*cP.charge()/constant::electromagnetic::epsilon0.value()/cP.mass());//Calculate the plasma frequency

    Info << "    Average plasma frequency        = " << plasmaFrequency << " Hz" << nl;
    if(warnTimeStep_ && plasmaFrequency*cloud.mesh().time().deltaTValue() > 2.0)
        Info << "    WARNING: Time step times plasma frequency is larger than 2!" << nl;

   //Reset list
   nParticle_ = 0.0;
}

// ************************************************************************* //
