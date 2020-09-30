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

#include "SimpleInjection.H"
#include "constants.H"
#include "triPointRef.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SimpleInjection<CloudType>::SimpleInjection
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud, typeName, associatedPatches),
    emitter_(cloud),
    injectorPatch_(-1),
    injectionSpeciesList_(this->coeffDict().lookup("injectionTypeIds")),
    injectionSpeciesTypeIds_(),
    interactionList_(),
    iF_(this->coeffDict().lookup("frequencies")),
    particleRemainder_(),
    eV_(this->coeffDict().lookup("eV")),
    reflectWithInitVel_(),
    isFloatingBC_(false),
    Qconv_(0.0),
    patchArea_(),
    boundaryCondition_(nullptr),
    injectionEndTime_(readScalar(this->coeffDict().lookup("injectionEndTime")))
{
    if(injectionSpeciesList_.size() != iF_.size() || iF_.size() != eV_.size())
    {
        FatalErrorInFunction
                << "not the same number of entries in frequencies, eV, injectionTypeIds" << nl
                << abort(FatalError);
    }
    particleRemainder_.setSize(injectionSpeciesList_.size(),0.0);

    interactionList_.setSize(cloud.typeIdList().size(),ptDeletionPatch);
    reflectWithInitVel_.setSize(cloud.typeIdList().size(),false);

    if(associatedPatches.size() > 1)
    {
        FatalErrorInFunction
                << "model does not support mutiple patches" << nl
                << abort(FatalError);
    }

    injectorPatch_ = associatedPatches[0];

    word velocityModel = this->coeffDict().lookup("velocityModel");
    if(velocityModel == "birdMaxwellianFlux")
        emitter_.initilizeParticleEmitter(injectorPatch_, ParticleEmitter<CloudType>::vmBirdMaxwellianFlux);
    else if(velocityModel == "mostProbableSpeed")
        emitter_.initilizeParticleEmitter(injectorPatch_, ParticleEmitter<CloudType>::vmMostProbableSpeed);
    else if(velocityModel == "maxwellianFlux")
        emitter_.initilizeParticleEmitter(injectorPatch_, ParticleEmitter<CloudType>::vmMaxwellianFlux);
    else if(velocityModel == "halfMaxwellianFlux")
        emitter_.initilizeParticleEmitter(injectorPatch_, ParticleEmitter<CloudType>::vmHalfMaxwellianFlux);
    else
    {
        FatalErrorInFunction
                << "valid injection velocity models are: " << nl
                << "birdMaxwellianFlux" << nl
                << "mostProbableSpeed" << nl
                << "maxwellianFlux" << nl
                << "halfMaxwellianFlux" << nl
                << abort(FatalError);
    }
    const polyPatch& patch = cloud.mesh().boundaryMesh()[injectorPatch_];

    if (!isType<polyPatch>(patch))
    {
        FatalErrorInFunction << patch.name() << " is not of type patch" << abort(FatalError);
    }

    const dictionary& interactions(this->coeffDict().subDict("interactions"));

     forAllConstIter(IDLList<entry>, interactions, iter)
     {
         if(iter().isDict())
         {
             const dictionary& subDict = iter().dict();
             label id = findIndex(this->owner().typeIdList(),iter().keyword());
             if(id == -1)
                 FatalErrorInFunction << "Species " << iter().keyword() << " is not defined in cloud" << endl;

             word interactionType = subDict.lookup("type");
             reflectWithInitVel_[id] = subDict.lookupOrDefault<bool>("useInitalVelocity",false);

             if(interactionType == "deletion")
                 interactionList_[id] = ptDeletionPatch;
             else if(interactionType == "reflection")
                 interactionList_[id] = ptReflection;
             else if(interactionType == "floating") {
                 interactionList_[id] = ptFloating;
                 isFloatingBC_ = true;
             }
             else
                 FatalErrorInFunction << "interaction type " << interactionType << " not defined. valid options [deletion/reflection/floating]" << abort(FatalError);
         }
     }

     if(isFloatingBC_)
     {
         volScalarField& phiE = cloud.elpotentialField();

         if(!isA<circuitBoundaryFvPatchField>(phiE.boundaryField()[injectorPatch_]))
             FatalErrorInFunction << "Expected circuitBoundary boundary condition for field " << phiE.name() << " on patch " << patch.name() << nl << abort(FatalError);

         boundaryCondition_ = &(refCast<circuitBoundaryFvPatchField>(phiE.boundaryFieldRef()[injectorPatch_]));

         const scalarField magSf(mag(patch.faceAreas()));
         patchArea_ = sum(magSf);
         reduce(patchArea_, sumOp<scalar>());
     }


     DynamicList<label> typeIds;
     forAll(injectionSpeciesList_, i)
     {
         word type = injectionSpeciesList_[i];
         label typeId = findIndex(cloud.typeIdList(), type);

         if (typeId == -1)
         {
             FatalErrorInFunction
                     << "species " << type << " not defined in cloud." << nl
                     << abort(FatalError);
         }
         typeIds.append(typeId);
     }
     injectionSpeciesTypeIds_.transfer(typeIds);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SimpleInjection<CloudType>::~SimpleInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::SimpleInjection<CloudType>::preUpdate_Boundary()
{
    if (isFloatingBC_)
    {
        //Circuit update
        reduce(Qconv_, sumOp<scalar>());

        boundaryCondition_->circuitGradient() = boundaryCondition_->circuitGradient() + Qconv_/patchArea_/constant::electromagnetic::epsilon0.value();

        Qconv_ = 0.0;
    }
}

template<class CloudType>
void Foam::SimpleInjection<CloudType>::injection()
{
      CloudType& cloud(this->owner());
      const fvMesh& mesh(cloud.mesh());


      if(mesh.time().value() > injectionEndTime_)
      {
          Info<< "[" << typeName << "] no injection (t > injectionEndTime_)" << endl;
          return;
      }
      const scalar dt = mesh.time().deltaTValue();

      forAll(injectionSpeciesTypeIds_, injS)
      {
          label injTypeId = injectionSpeciesTypeIds_[injS];
          label particlesInserted = 0;

          scalar nParcels =  iF_[injS]/cloud.constProps(injTypeId).nParticle() * dt + particleRemainder_[injS];//if_ = frequency/nParticle

          label nParcelsToInject = label(nParcels);

          particleRemainder_[injS] = nParcels-nParcelsToInject;

          for(label n = 0; n < nParcelsToInject; n++)
          {
              scalar temperature = eV_[injS]*constant::electromagnetic::e.value()/constant::physicoChemical::k.value();
              typename CloudType::parcelType* parcel = emitter_.emitParticle(temperature,injTypeId);
              if(parcel != nullptr)
                  particlesInserted++;
          }

          reduce(particlesInserted, sumOp<label>());

          Info<< "[" << typeName << "] Particles inserted (" << injectionSpeciesList_[injS] << ")           = "
                  << particlesInserted << " | " << particleRemainder_[injS] << endl;
      }
}


template<class CloudType>
bool Foam::SimpleInjection<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    label typeId = p.typeId();
    if(injectorPatch_ == p.patch())
    {
        if(interactionList_[typeId] == ptReflection)
        {
            if(!reflectWithInitVel_[typeId]) {
                p.specularReflect();//do not update fields
            }
            else
            {


                label injS = findIndex(injectionSpeciesTypeIds_,typeId);

                if(injS != -1)
                {
                    scalar temperature = eV_[injS]*constant::electromagnetic::e.value()/constant::physicoChemical::k.value();
                    td.keepParticle = false;
                    typename CloudType::parcelType* newp = emitter_.emitParticleExplicit(temperature, typeId, p.stepFraction());
                    newp->nParticle() = p.nParticle();
                    //FIXME: this should also set the chargeModifier
                }
                else//not injected so normal reflection
                {
                    p.specularReflect();//do not update fields
                }
            }
            return true;
        }
        else if(interactionList_[typeId] == ptDeletionPatch)
        {
            td.keepParticle = false;
            return true;
        }
        else if(interactionList_[typeId] == ptFloating)
        {
            const scalar charge = p.nParticle()*p.charge();
            Qconv_ += charge;

            td.keepParticle = false;//simply delete particle do not update fields
            return true;
        }
    }


    return false;
}

template<class CloudType>
void Foam::SimpleInjection<CloudType>::particleEjection(typename CloudType::parcelType& p, label patchId)
{
    if(isFloatingBC_){
        scalar Q = p.charge()*p.nParticle();
        Qconv_ -= Q;
    }
}

// ************************************************************************* //
