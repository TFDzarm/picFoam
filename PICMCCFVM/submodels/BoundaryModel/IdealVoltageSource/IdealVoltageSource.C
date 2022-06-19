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

#include "IdealVoltageSource.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"
#include "fixedValueFvPatchField.H"
#include "processorPolyPatch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IdealVoltageSource<CloudType>::IdealVoltageSource
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud, typeName,associatedPatches),
    interactionList_(),
    chargeAccumulator_(0.0),
    chargeAccumulatorCathode_(0.0),
    anodeAverage_(0.0),
    cathodeAverage_(0.0),
    averageStart_(cloud.mesh().time().value()),
    Qconv_(0.0),
    Qconv1_(0.0),
    Qconv2_(0.0),
    QconvCathode_(0.0),
    Qconv1Cathode_(0.0),
    Qconv2Cathode_(0.0),
    emissionList_(this->coeffDict(),cloud)
{
    //Get number of patches not counting processorPatches
    label nPatches = cloud.mesh().boundaryMesh().size();
    if(Pstream::parRun())
    {
        forAll(cloud.mesh().boundaryMesh(), patchId)
        {
            const polyPatch& patch = cloud.mesh().boundaryMesh()[patchId];
            if (isA<processorPolyPatch>(patch))
            {
                nPatches--;
            }
        }
    }

    interactionList_.setSize(nPatches,etNone);

    volScalarField& phiE = cloud.elpotentialField();
    //Check patches
    forAll(associatedPatches, i)
    {
        label patchId = associatedPatches[i];
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patchId];

        if(isA<processorPolyPatch>(patch))
            continue;//Could we break?

        if (isType<polyPatch>(patch))
        {
            const dictionary& patchPropertiesDict = this->coeffDict().subDict(patch.name());
            const word type(patchPropertiesDict.lookup("type"));

            if(type == "anode")
            {
                interactionList_[patchId] = etAnode;

                //Has to be of type fixedValue
                if(!isA<fixedValueFvPatchField<scalar>>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected fixedValue boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);
            }
            else if(type == "cathode")
            {
                emissionList_.initilizeAll(patchId,ParticleEmitter<CloudType>::vmMaxwellianFlux);
                interactionList_[patchId] = etCathode;

                //Has to be of type fixedValue
                if(!isA<fixedValueFvPatchField<scalar>>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected fixedValue boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);
            }
            else
            {
                FatalErrorInFunction
                        << "patch type " << type << " for patch " << patch.name() << " not defined. valid options: [anode/cathode]" << nl
                        << abort(FatalError);
            }
        }
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IdealVoltageSource<CloudType>::~IdealVoltageSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::IdealVoltageSource<CloudType>::postUpdate_Boundary()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    scalar dt = mesh.time().deltaTValue();

    reduce(chargeAccumulator_, sumOp<scalar>());
    reduce(chargeAccumulatorCathode_, sumOp<scalar>());

    //Charge that was accumulated previous timesteps
    Qconv2_ = Qconv1_;
    Qconv1_ = Qconv_;
    Qconv_ -= (chargeAccumulator_);

    Qconv2Cathode_ = Qconv1Cathode_;
    Qconv1Cathode_ = QconvCathode_;
    QconvCathode_ -= (chargeAccumulatorCathode_);

    chargeAccumulator_ = 0.0;
    chargeAccumulatorCathode_ = 0.0;

    //Calculate the running average of the current flowing through the circuit
    scalar beta = dt/(mesh.time().value()-averageStart_);

    anodeAverage_ = (1.0-beta)*anodeAverage_ + beta*((2.0*(Qconv_-Qconv1_) + 0.5*(Qconv2_-Qconv_))/dt);
    cathodeAverage_ = (1.0-beta)*cathodeAverage_ + beta*((2.0*(QconvCathode_-Qconv1Cathode_) + 0.5*(Qconv2Cathode_-QconvCathode_))/dt);

    Info << "[" << typeName << "] anode current: " << anodeAverage_ << " A" << nl
         << "[" << typeName << "] cathode current: " << cathodeAverage_ << " A" << endl;
}

/*
Foam::IdealVoltageSource<CloudType>::injection called in evolve function
*/
template<class CloudType>
void Foam::IdealVoltageSource<CloudType>::injection()
{
    //Tell EmissionModels to inject particles
    emissionList_.emission();
}


template<class CloudType>
void Foam::IdealVoltageSource<CloudType>::particleEjection(typename CloudType::parcelType& p, label patchId)
{
    //Have we ejected a charged particle? Subtract the charge...
    if(interactionList_[patchId] == etAnode)
        chargeAccumulator_ -= p.charge()*p.nParticle();
    else if(interactionList_[patchId] == etCathode)
        chargeAccumulatorCathode_ -= p.charge()*p.nParticle();
}

template<class CloudType>
bool Foam::IdealVoltageSource<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    CloudType& cloud(this->owner());

    label patchId = p.patch();
    const scalar charge = p.nParticle()*p.charge();

    //Did a charged particle hit the boundary? Add the charge... Else reflect...
    if(charge != 0.0)
    {
        if(interactionList_[patchId] == etAnode)
            chargeAccumulator_ += charge;
        else if(interactionList_[patchId] == etCathode)
            chargeAccumulatorCathode_ += charge;

        p.wallAbsorption(cloud, td);
    }
    else {
        p.wallReflection(cloud, td);
    }

    //Check EmissionModels
    if(interactionList_[patchId] == etCathode)
        emissionList_.collisionalEmission(p,td);

    return true;
}

// ************************************************************************* //
