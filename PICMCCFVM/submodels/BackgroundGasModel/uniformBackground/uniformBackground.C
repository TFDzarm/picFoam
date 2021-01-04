/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "uniformBackground.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
Foam::uniformBackground<CloudType>::uniformBackground

Read in the options from constant/picProperties
To define the density one can specify the number Density directly or the pressure.
*/
template<class CloudType>
Foam::uniformBackground<CloudType>::uniformBackground
(
    const dictionary& dict,
    CloudType& owner
)
:
    BackgroundGasModel<CloudType>(dict,owner,typeName),
    species_(-1),
    numberDensity_(owner.mesh().cells().size(),0.0),
    temperature_(owner.mesh().cells().size(),0.0),
    nParticle_(readScalar(this->coeffDict().lookup("nEquivalentParticles"))),
    velocity_(this->coeffDict().lookup("velocity"))//drift velocity
{

    word species = this->coeffDict().lookup("species");
    label id = findIndex(owner.typeIdList(),species);
    if(id < 0)
        FatalErrorInFunction << "Undefined typeId " << species << abort(FatalError);

    scalar numberDensity = this->coeffDict().lookupOrDefault("numberDensity",-1.0);//Number density, default value if entry doesn't exist -1
    scalar temperature = readScalar(this->coeffDict().lookup("temperature"));
    temperature_ = temperature;

    species_ = id;
    if(numberDensity < 0.0)//If the default value is set, check if pressure was given
    {
        scalar pressure = this->coeffDict().lookupOrDefault("pressure",-1.0);
        if(pressure < 0.0)
            FatalErrorInFunction << "No keyword pressure or numberDensity defined." << abort(FatalError);

        if(temperature > VSMALL)
            numberDensity_ = pressure/(temperature_*constant::physicoChemical::k.value());//assume ideal gas law
        else
            numberDensity_ = pressure/(VSMALL*constant::physicoChemical::k.value());
    }
    else
        numberDensity_ = numberDensity;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::uniformBackground<CloudType>::~uniformBackground()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::uniformBackground<CloudType>::active() const
{
    return true;
}

template<class CloudType>
Foam::label Foam::uniformBackground<CloudType>::species() const
{
    return species_;
}



template<class CloudType>
Foam::scalar Foam::uniformBackground<CloudType>::nParticle() const
{
    return nParticle_;
}

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::uniformBackground<CloudType>::numberDensity() const
{
    tmp<scalarField> n(numberDensity_);//Return a tmp field
    return n;
}

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::uniformBackground<CloudType>::temperature() const
{
    tmp<scalarField> T(temperature_);//Return a tmp field
    return T;
}

template<class CloudType>
Foam::vector Foam::uniformBackground<CloudType>::sampleVelocity(label celli)
{
    CloudType& cloud(this->owner());

    scalar mass = cloud.constProps(species_).mass();
    return cloud.equipartitionLinearVelocity(temperature_[celli],mass)+velocity_;//Sample a velocity from Maxwell-Boltzmann distribution and add the drift velocity
}


// ************************************************************************* //
