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

#include "circuitRLC.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"
#include "fixedGradientFvPatchField.H"
#include "fixedValueFvPatchField.H"
#include "processorPolyPatch.H"
#include "EmissionModel.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::circuitRLC<CloudType>::circuitRLC
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud, typeName,associatedPatches),
    interactionList_(),
    sigma1_(0.0),
    chargeAccumulator_(0.0),
    cathodePatch_(-1),
    anodePatch_(-1),
    cathodeArea_(),
    boundaryCondition_(nullptr),
    V_(readScalar(this->coeffDict().lookup("V"))),
    omega_(this->coeffDict().lookupOrDefault("omega",0.0)),
    phase_(this->coeffDict().lookupOrDefault("phase",constant::mathematical::piByTwo)),
    Qex_(0.0),//this->coeffDict().lookupOrDefault("Q0",0.0)
    Qex1_(0.0),
    Qex2_(0.0),
    Qex3_(0.0),
    R_(readScalar(this->coeffDict().lookup("R"))),
    C_(readScalar(this->coeffDict().lookup("C"))),
    L_(readScalar(this->coeffDict().lookup("L"))),
    a0_(0.0),
    a1_(0.0),
    a2_(0.0),
    a3_(0.0),
    a4_(0.0),
    k_(0.0),
    Iex_(0.0),
    epsilon_(this->coeffDict().lookupOrDefault("epsilon",constant::electromagnetic::epsilon0.value())),
    emissionList_(this->coeffDict(),cloud)
{
    scalar dt = cloud.mesh().time().deltaTValue();
    a0_ = 2.25*L_/dt/dt + 1.5*R_/dt + 1.0/C_;
    a1_ = -6.0*L_/dt/dt - 2.0*R_/dt;
    a2_ = 5.5*L_/dt/dt + 0.5*R_/dt;
    a3_ = -2.0*L_/dt/dt;
    a4_ = 0.25*L_/dt/dt;

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

                if(!isA<fixedValueFvPatchField<scalar>>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected fixedValue boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

            }
            else if(type == "cathode")
            {
                emissionList_.initilizeAll(patchId,ParticleEmitter<CloudType>::vmMaxwellianFlux);

                interactionList_[patchId] = etCathode;

                cathodePatch_ = patchId;


                if(!isA<RLCBoundaryFvPatchField>(phiE.boundaryField()[patchId]))
                    FatalErrorInFunction << "Expected RLCBoundary boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

                const scalarField magSf(mag(patch.faceAreas()));
                cathodeArea_ = sum(magSf);
                reduce(cathodeArea_, sumOp<scalar>());

                boundaryCondition_ = &(refCast<RLCBoundaryFvPatchField>(phiE.boundaryFieldRef()[patchId]));
                Qex_ = boundaryCondition_->Q();
                Qex1_ = boundaryCondition_->Q1();
                Qex2_ = boundaryCondition_->Q2();
                Qex3_ = boundaryCondition_->Q3();
                sigma1_ = boundaryCondition_->sigma();
                boundaryCondition_->coeff1() = 1.0/(a0_*epsilon_*cathodeArea_+1.0/boundaryCondition_->patch().deltaCoeffs());
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
Foam::circuitRLC<CloudType>::~circuitRLC()
{


}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::circuitRLC<CloudType>::preUpdate_Boundary()
{
    reduce(chargeAccumulator_, sumOp<scalar>());

    k_ = a1_*Qex_+a2_*Qex1_+a3_*Qex2_+a4_*Qex3_;

    boundaryCondition_->coeff2() = a0_*cathodeArea_*sigma1_+a0_*(chargeAccumulator_-Qex_)+(source()-k_);
}

template<class CloudType>
void Foam::circuitRLC<CloudType>::postUpdate_Boundary()
{

    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    scalar dt = mesh.time().deltaTValue();

    k_ = (a1_*Qex_ + a2_*Qex1_ + a3_*Qex2_ + a4_*Qex3_)/a0_;

    scalar phiBound = 0.0;

    //Calc an average of the boundary potential (there is probably a better way of doing this)
    if(Pstream::parRun())
    {
        label count = 0;
        phiBound = sum(*boundaryCondition_);
        count = (*boundaryCondition_).size();

        reduce(phiBound, sumOp<scalar>());
        reduce(count, sumOp<label>());
        phiBound /= count;

    }
    else
        phiBound = voltage();

    Qex3_ = Qex2_;
    Qex2_ = Qex1_;
    Qex1_ = Qex_;
    Qex_ = (source()-phiBound)/a0_-k_;

    Iex_ = (Qex_-Qex1_)/dt;
    sigma1_ += (chargeAccumulator_+Iex_*dt)/cathodeArea_;

    boundaryCondition_->Q() = Qex_;
    boundaryCondition_->Q1() = Qex1_;
    boundaryCondition_->Q2() = Qex2_;
    boundaryCondition_->Q3() = Qex3_;
    boundaryCondition_->sigma() = sigma1_;

    chargeAccumulator_ = 0.0;

    //sanity check
    //scalar check = L_*(2.0*Qex_-5.0*Qex1_+4.0*Qex2_-Qex3_)/dt/dt+R_*(3.0*Qex_-4.0*Qex1_+Qex2_)/(2.0*dt);
    //Info << "circ: " << check << " == " << V_-(*boundaryCondition_)[0] << endl;

    Info << "[" << typeName << "] voltage: " << phiBound << " V current: " << Iex_ << " A" << endl;
}

template<class CloudType>
void Foam::circuitRLC<CloudType>::injection()
{
    emissionList_.emission();
}


template<class CloudType>
void Foam::circuitRLC<CloudType>::particleEjection(typename CloudType::parcelType& p, label patchId)
{
    if(interactionList_[patchId] == etCathode)
    {
        scalar Q = p.charge()*p.nParticle();
        chargeAccumulator_ -= Q;
    }
}

template<class CloudType>
bool Foam::circuitRLC<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    CloudType& cloud(this->owner());

    label patchId = p.patch();
    const scalar charge = p.nParticle()*p.charge();

    if(charge != 0.0)
    {
        if(interactionList_[patchId] == etCathode)
            chargeAccumulator_ += charge;
        p.wallAbsorption(cloud, td);
    }
    else {
        p.wallReflection(cloud, td);
    }

    if(interactionList_[patchId] == etCathode)
        emissionList_.collisionalEmission(p,td);

    return true;
}

// ************************************************************************* //
