/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "IdealCurrentSource.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"
#include "fixedGradientFvPatchField.H"
#include "fixedValueFvPatchField.H"
#include "processorPolyPatch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IdealCurrentSource<CloudType>::IdealCurrentSource
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud, typeName,associatedPatches),
    interactionList_(),
    sigma_(0.0),
    chargeAccumulator_(0.0),
    cathodePatch_(-1),
    anodePatch_(-1),
    anodeArea_(),
    cathodeArea_(),
    boundaryCondition_(nullptr),
    I_(this->coeffDict().lookupOrDefault("I",-1.0)),
    J_(this->coeffDict().lookupOrDefault("J",-1.0)),
    emissionList_(this->coeffDict(),cloud)
{
    if(I_ < 0.0 && J_ < 0.0)
    {
        FatalErrorInFunction << "No current \"I\" or current density \"J\" was given." << abort(FatalError);
    }

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
                anodePatch_ = patchId;

                //The anode has to be a fixedValue boundary condition
                if(!isA<fixedValueFvPatchField<scalar>>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected fixedValue boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

                //Calculate area if the patch
                const scalarField magSf(mag(patch.faceAreas()));
                anodeArea_ = sum(magSf);
                reduce(anodeArea_, sumOp<scalar>());

            }
            else if(type == "cathode")
            {
                emissionList_.initilizeAll(patchId,ParticleEmitter<CloudType>::vmMaxwellianFlux);//Allways use a MaxwellianFlux emitter

                interactionList_[patchId] = etCathode;

                cathodePatch_ = patchId;

                //Cathode has to be a circuit boundary
                if(!isA<circuitBoundaryFvPatchField>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected circuitBoundary boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

                boundaryCondition_ = &(refCast<circuitBoundaryFvPatchField>(phiE.boundaryFieldRef()[patchId]));

                //Calculate area if the patch
                const scalarField magSf(mag(patch.faceAreas()));
                cathodeArea_ = sum(magSf);
                reduce(cathodeArea_, sumOp<scalar>());

                //If the real area is larger e.g. we run a reduced mesh case -> scale the current
                scalar refArea = this->coeffDict().lookupOrDefault("refArea",0.0);
                if(refArea > 0.0) {
                    I_ *= cathodeArea_/refArea;
                    Info << "Cathode area: " << cathodeArea_ << " refArea: " << refArea << " setting I to " << I_ << endl;
                }
                //If the current density J instead of I was give set I
                if(J_ > 0.0)
                {
                    I_ = J_*cathodeArea_;
                    Info << "Current density J=" << J_ << " given, set I to " << I_ << endl;
                }

                //Initial surface charge
                sigma_ = boundaryCondition_->circuitGradient()*constant::electromagnetic::epsilon0.value();

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
Foam::IdealCurrentSource<CloudType>::~IdealCurrentSource()
{


}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IdealCurrentSource<CloudType>::preUpdate_Boundary()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    //Parallel COM
    reduce(chargeAccumulator_, sumOp<scalar>());

    scalar dt = mesh.time().deltaTValue();

    //Update the surface charge
    sigma_ += (chargeAccumulator_ - I_*dt)/cathodeArea_;

    //Update the gradient in the boundary condition
    boundaryCondition_->circuitGradient() = sigma_/constant::electromagnetic::epsilon0.value();
    Info << "[" << typeName << "] surface charge: " << sigma_*cathodeArea_  << " As" << endl;

    //Clear the charge accumulator
    chargeAccumulator_ = 0.0;
}

template<class CloudType>
void Foam::IdealCurrentSource<CloudType>::postUpdate_Boundary()
{
    //Calculate the average voltage on the boundary
    scalar avgV = 0.0;
    if(boundaryCondition_->size() > 0) { //In parallel runs some processors might not have a part of the boundary
        const scalarField& faceAreas = boundaryCondition_->patch().magSf();
        avgV = sum(faceAreas*boundaryCondition_->patchInternalField());
    }
    reduce(avgV, sumOp<scalar>());
    Info << "[" << typeName << "] average voltage: " << avgV/cathodeArea_ << " V " << endl;
}

/*
Foam::IdealCurrentSource<CloudType>::injection called in evolve function
*/
template<class CloudType>
void Foam::IdealCurrentSource<CloudType>::injection()
{
    //Tell EmissionModels to inject particles
    emissionList_.emission();
}


template<class CloudType>
void Foam::IdealCurrentSource<CloudType>::particleEjection(typename CloudType::parcelType& p, label patchId)
{
    //Have we ejected a charged particle? Subtract the charge...
    if(interactionList_[patchId] == etCathode)
    {
        scalar Q = p.charge()*p.nParticle();
        chargeAccumulator_ -= Q;
    }
}

template<class CloudType>
bool Foam::IdealCurrentSource<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    CloudType& cloud(this->owner());

    label patchId = p.patch();
    const scalar charge = p.nParticle()*p.charge();

    //Did a charged particle hit the boundary? Add the charge... Else reflect...
    if(charge != 0.0)
    {
        if(interactionList_[patchId] == etCathode)
            chargeAccumulator_ += charge;

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
