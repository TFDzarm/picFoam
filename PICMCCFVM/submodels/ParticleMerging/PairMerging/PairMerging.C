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

#include "fvCFD.H"
#include "PairMerging.H"
#include "Random.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairMerging<CloudType>::PairMerging
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ParticleMerging<CloudType>(dict,cloud,typeName),
    reductionFactor_(readScalar(this->coeffDict().lookup("reductionFactor"))),
    cellParcelMax_(readLabel(this->coeffDict().lookup("cellParcelsMax"))),
    maxChecks_(readLabel(this->coeffDict().lookup("maxChecks"))),
    nCheckNeighbors_(readLabel(this->coeffDict().lookup("nCheckNeighbors"))),
    angleDifference_(readScalar(this->coeffDict().lookup("angleDifference"))*constant::mathematical::pi/180.0),
    velocityRatio_(readScalar(this->coeffDict().lookup("velocityRatio"))),
    cumulativeKinEnergyError_(0.0),
    mergeSpecies_(),
    nParticleDev_(this->coeffDict().lookup("nDeviation"))
{
    if((cellParcelMax_*reductionFactor_) < nCheckNeighbors_)
        FatalErrorInFunction << "Merging: Reduction to " << cellParcelMax_*reductionFactor_ << " while checking " << nCheckNeighbors_ << " neighbors is not possible" << abort(FatalError);

    DynamicList<label> mS;
    wordList species(this->coeffDict().lookup("species"));
    forAll(species, i)
    {
        word type = species[i];
        label typeId = findIndex(cloud.typeIdList(), type);

        if (typeId == -1)
        {
            FatalErrorInFunction
                    << "species " << type << " not defined in cloud." << nl
                    << abort(FatalError);
        }
        mS.append(typeId);
    }
    mergeSpecies_.transfer(mS);

    Warning << "Particle merging is untested. Use with caution!" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairMerging<CloudType>::~PairMerging()
{}

template<class CloudType>
bool Foam::PairMerging<CloudType>::shouldMerge() const
{
    return true;//check every timestep
}

template<class CloudType>
void Foam::PairMerging<CloudType>::mergeParticles()
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

        //does not include newly created ions and electrons of this timestep!
        List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

        label mergedParcel = 0;
        label checkedParcel = 0;
        scalar kinEnergyError = 0;
        forAll(sortedCellOccupancy, celli)// Go through all cells
        {
            forAll(mergeSpecies_,idx)
            {
                label typeId = mergeSpecies_[idx];

                label N = sortedCellOccupancy[celli][typeId].size();

                if(N > cellParcelMax_)
                {
                    label nTarget = N*reductionFactor_;

                    List<scalar> mergeVelocities(N);
                    forAll(sortedCellOccupancy[celli][typeId],i)
                            mergeVelocities[i] = mag(sortedCellOccupancy[celli][typeId][i]->U());


                    List<label> mergeSelect(N);
                    forAll(mergeSelect,pos)//iota
                        mergeSelect[pos] = pos;

                    sort(mergeSelect,[&mergeVelocities](label i1, label i2) {return mergeVelocities[i1] < mergeVelocities[i2];});//Lambda magic


                    label checks = 0;
                    while(N > nTarget)
                    {
                        label pId1 = rndGen.sampleAB<label>(0, N);
                        typename CloudType::parcelType* p1 = sortedCellOccupancy[celli][typeId][mergeSelect[pId1]];
                        typename CloudType::parcelType* p2 = nullptr;

                        label pId2 = pId1-(nCheckNeighbors_/2);
                        if(pId2 < 0)
                            pId2 = 0;
                        if(pId2+nCheckNeighbors_ > N-1)
                            pId2 = N-nCheckNeighbors_-1;

                        for(label i = 0; i < nCheckNeighbors_; i++)
                        {
                            if(pId2==pId1)
                                pId2++;

                            typename CloudType::parcelType* candidate = sortedCellOccupancy[celli][typeId][mergeSelect[pId2]];

                            //check
                            checkedParcel++;

                            if(nParticleDev_[idx] - mag(candidate->nParticle() - p1->nParticle()) < 0.0) {
                                pId2++;
                                continue;
                            }

                            if(candidate->chargeModifier() != p1->chargeModifier()) {
                                pId2++;
                                continue;
                            }


                            if((candidate->U()&p1->U()) > mergeVelocities[pId2]*mergeVelocities[pId1]*::cos(angleDifference_)) {
                                pId2++;
                                continue;
                            }

                            scalar ratio = mergeVelocities[pId2]/mergeVelocities[pId1];
                            if(ratio > 1.0)
                                ratio = 1.0/ratio;
                            if(ratio < velocityRatio_) {
                                pId2++;
                                continue;
                            }

                            p2 = candidate;
                            break;

                        }

                        if(p2 != nullptr)
                        {
                            mergedParcel++;
                            scalar nSumParticle = p1->nParticle()+p2->nParticle();
                            const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(typeId);

                            vector CoM = (p1->nParticle()*p1->position() + p2->nParticle()*p2->position())/nSumParticle;
                            //FIXME: account for displacement in potential field?
                            vector newU = (p1->nParticle()*p1->U() + p2->nParticle()*p2->U())/nSumParticle;

                            //Error checking
                            scalar preKinE = 0.5*cP.mass()*( (p1->U() & p1->U())*p1->nParticle() + (p2->U() & p2->U())*p2->nParticle());
                            scalar postKinE = 0.5*cP.mass()*(newU&newU)*nSumParticle;
                            kinEnergyError += postKinE-preKinE;


                            //Apply

                            //Does not work
                            //p2->position() = CoM; // this does not do anything!
                            //p2->locate(CoM,nullptr,p2->cell(),false,""); // is private
                            //p2->U() = newU;

                            //Create a new merged parcel instead
                            typename CloudType::parcelType* p = cloud.addNewParcel(CoM,p2->cell(),newU,p2->typeId());
                            p->nParticle() = nSumParticle;

                            cloud.deleteParticle(*p2);
                            sortedCellOccupancy[celli][typeId][mergeSelect[pId2]] = nullptr;
                            cloud.deleteParticle(*p1);
                            sortedCellOccupancy[celli][typeId][mergeSelect[pId1]] = nullptr;

                            label minPos,maxPos;
                            if(pId1 > pId2) {
                                minPos = pId2;maxPos = pId1;
                            } else {
                                minPos = pId1;maxPos = pId2;
                            }

                            //shift left
                            for(label i = minPos; i < N-2; i++)
                            {
                                if(i >= maxPos-1)
                                    mergeSelect[i] = mergeSelect[i+2];
                                else
                                    mergeSelect[i] = mergeSelect[i+1];
                            }
                            N-=2;
                        }
                        checks++;
                        if(checks > maxChecks_)
                            break;
                    }
                }
            }
        }
        reduce(mergedParcel,sumOp<label>());
        reduce(checkedParcel,sumOp<label>());

        reduce(kinEnergyError,sumOp<scalar>());

        cumulativeKinEnergyError_ += kinEnergyError;

        if(mergedParcel > 0)
        {
            Info << "    Merged parcels                          = " << mergedParcel << nl
                 << "    Acceptance rate                         = " << scalar(mergedParcel)/scalar(checkedParcel) << nl
                 << "    Kinetic energy error (cumulative error) = " << kinEnergyError << " (" << cumulativeKinEnergyError_ << ")" << endl;
        }
        else
            Info << "    No merging" << endl;
}

// ************************************************************************* //
