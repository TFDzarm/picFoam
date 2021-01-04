/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "BinaryCollisionModel.H"
#include "TotalCrossSectionModel.H"
#include "BackgroundGasModel.H"
#include "WeightCorrectionModel.H"
#include "HardSphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BinaryCollisionModel<CloudType>::BinaryCollisionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    sigmaTcRMax_//Read the field from the time dir
    (
        IOobject
        (
            "picSigmaTcRMax",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.mesh()
    ),
    sigmaTIoncRMax_//Read the field from the time dir
    (
        IOobject
        (
            "picSigmaTcRMax_IonNeutral",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.mesh()
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            "pic:collisionSelectionRemainder",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("collisionSelectionRemainder", dimless, 0)
    ),
    collisionSelectionIonRemainder_
    (
        IOobject
        (
            "pic:collisionSelectionIonRemainder",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("collisionSelectionIonRemainder", dimless, 0)
    ),
    totalCrossSection_(),
    backgroundCollisionPropIon_
    (
        IOobject
        (
            "backgroundCollisionProp_IonNeutral",
             owner.mesh().time().timeName(),
             owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionProp_IonNeutral", dimless, 0)
    ),
    backgroundCSRIon_
    (
        IOobject
        (
            "backgroundCollisionRemainder_IonnNeutral",
             owner.mesh().time().timeName(),
             owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_IonNeutral", dimless, 0)
   ),
   backgroundCollisionProp_
   (
       IOobject
       (
           "backgroundCollisionProp_NeutralNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
       ),
       owner.mesh(),
       dimensionedScalar("backgroundCollisionProp_NeutralNeutral", dimless, 0)
   ),
   backgroundCSR_
   (
       IOobject
       (
           "backgroundCollisionRemainder_NeutralNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
       ),
       owner.mesh(),
       dimensionedScalar("backgroundCollisionRemainder_NeutralNeutral", dimless, 0)
   ),
   weightCorrection_(),
   chargeExchangeCollision_(false)
{}


template<class CloudType>
Foam::BinaryCollisionModel<CloudType>::BinaryCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    sigmaTcRMax_//Read the field from the time dir
    (
        IOobject
        (
            "picSigmaTcRMax",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.mesh()
    ),
    sigmaTIoncRMax_//Read the field from the time dir
    (
        IOobject
        (
            "picSigmaTcRMax_IonNeutral",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.mesh()
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            "pic:collisionSelectionRemainder",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("collisionSelectionRemainder", dimless, 0)
    ),
    collisionSelectionIonRemainder_
    (
        IOobject
        (
            "pic:collisionSelectionIonRemainder",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("collisionSelectionIonRemainder", dimless, 0)
    ),
    totalCrossSection_//Construct the model
    (
        TotalCrossSectionModel<CloudType>::New
        (
            coeffDict_,
            owner
        )
    ),
    backgroundCollisionPropIon_
    (
        IOobject
        (
            "backgroundCollisionProp_IonNeutral",
             owner.mesh().time().timeName(),
             owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionProp_IonNeutral", dimless, 0)
    ),
    backgroundCSRIon_
    (
        IOobject
        (
            "backgroundCollisionRemainder_IonnNeutral",
             owner.mesh().time().timeName(),
             owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_IonNeutral", dimless, 0)
   ),
   backgroundCollisionProp_
   (
       IOobject
       (
           "backgroundCollisionProp_NeutralNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
       ),
       owner.mesh(),
       dimensionedScalar("backgroundCollisionProp_NeutralNeutral", dimless, 0)
   ),
   backgroundCSR_
   (
       IOobject
       (
           "backgroundCollisionRemainder_NeutralNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
       ),
       owner.mesh(),
       dimensionedScalar("backgroundCollisionRemainder_NeutralNeutral", dimless, 0)
   ),
   weightCorrection_(//Construct the model
       WeightCorrectionModel<CloudType>::New
       (
           coeffDict_,
           owner
       )
   ),
   chargeExchangeCollision_(readBool(this->coeffDict().lookup("handleChargeExchange")))//Do we use the simple model descriped by Nanbu?
{
    if(chargeExchangeCollision_ && !isA<HardSphere<CloudType>>(totalCrossSection_()))
    {
        FatalErrorInFunction << "Charge exchange can only be used in combination with the HardSphere cross section model." << abort(FatalError);
    }
    if(chargeExchangeCollision_ && owner.neutralSpecies().size() > 1)
    {
        FatalErrorInFunction << "Charge exchange is only valid for one neutral and its ion species." << abort(FatalError);
    }

    const BackgroundGasModel<CloudType>& backgroundGas(owner.backgroundGas());
    if(backgroundGas.active())//If we use a background gas model calculate the collision probability
    {
        label bgTypeId = backgroundGas.species();
        const scalarField& bgNumberDensity = backgroundGas.numberDensity();

        scalar dt = owner.mesh().time().deltaTValue();

        //Use sigmaTIoncRMax_ as inital value
        forAll(backgroundCollisionPropIon_,i)
            backgroundCollisionPropIon_[i] = 1.0 -::exp(-sigmaTIoncRMax_[i]*bgNumberDensity[i]*dt);

        //Use sigmaTcRMax_ as inital value
        forAll(backgroundCollisionProp_,i)
            backgroundCollisionProp_[i] = 1.0 -::exp(-sigmaTcRMax_[i]*bgNumberDensity[i]*dt);

        Info << "|->    BackgroundGasModel " << backgroundGas.type() << " neutral collision properties:" << nl
             << "|          species: " << owner.typeIdList()[bgTypeId] << nl
             << "|          avg number density: " << average(bgNumberDensity) << " m^-3" << nl
             << "|          avg temperature: " << average(backgroundGas.temperature()) << " K" << nl
             << "|          avg sigmaTcRmax: " << average(sigmaTcRMax_.primitiveField()) << nl
             << "|          avg collision Pmax: " << average(backgroundCollisionProp_.field()) << nl

             << "|->    BackgroundGasModel " << backgroundGas.type() << " ion collision properties:" << nl
             << "           avg sigmaTcRmax: " << average(sigmaTIoncRMax_.primitiveField()) << nl
             << "           avg collision Pmax: " << average(backgroundCollisionPropIon_.field()) << nl << endl;
    }
    forAll(collisionSelectionRemainder_, i)//Initialize the remainders to a random value
    {
        collisionSelectionRemainder_[i] = owner.rndGen().scalar01();
        collisionSelectionIonRemainder_[i] = owner.rndGen().scalar01();
        backgroundCSRIon_[i] = owner.rndGen().scalar01();
        backgroundCSR_[i] = owner.rndGen().scalar01();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BinaryCollisionModel<CloudType>::~BinaryCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::BinaryCollisionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType&
Foam::BinaryCollisionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary&
Foam::BinaryCollisionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::BinaryCollisionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}

/*
Foam::BinaryCollisionModel<CloudType>::chargeExchangeCollision

Perform a charge exchange, use a simple correction method of Particle is different.
*/
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::chargeExchangeCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ)
{
    scalar rndU = this->owner().rndGen().scalar01();

    //Identity Switch
    vector tmpP = pP->U();
    if(rndU <= pQ->nParticle()/pP->nParticle())
        pP->U() = pQ->U();

    if(rndU <= pP->nParticle()/pQ->nParticle())
        pQ->U() = tmpP;
}

/*
Foam::BinaryCollisionModel<CloudType>::chargeExchangeCollision

Charge exchange collision with the background
*/
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::chargeExchangeCollision(typename CloudType::parcelType* pP, vector& Uq, label idQ)
{
    vector tmpP = pP->U();
    pP->U() = Uq;
    Uq = tmpP;
}

/*
Foam::BinaryCollisionModel<CloudType>::elasticCollision

Elastic collision with the background. Update the velocity by calling the chosen model.
*/
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::elasticCollision(typename CloudType::parcelType* pP, vector& Uq, label idQ)
{
    updateVelocity(*pP,Uq,idQ);
}

/*
Foam::BinaryCollisionModel<CloudType>::elasticCollision

Elastic collision call the chosen model and correct the velocity if the particle weights are different.
*/
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::elasticCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ)
{
    vector preUp = pP->U();
    vector preUq = pQ->U();

    updateVelocity(*pP,*pQ);
    weightCorrection_->correctVelocity(pP,pQ,preUp,preUq);
}


template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::performCollisions()
{

    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    const List<label>& chargedSpecies(cloud.chargedSpecies());
    scalar deltaT = mesh.time().deltaTValue();

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    //Used for counting the collision events
    label collisionCandidates = 0;
    label collisions = 0;

    label ionCollisionCandidates = 0;
    label ionCollisions = 0;

    label chargeExchanges = 0;

    DynamicList<typename CloudType::parcelType*> collisionList;//Lists of particles that can collide
    DynamicList<typename CloudType::parcelType*> ionCollisionList;

    forAll(sortedCellOccupancy, celli)// Go through all cells, add particles and look for the highest nParticle
    {
        label nC(0), iC(0);
        collisionList.clear();
        ionCollisionList.clear();
        scalar nParticleNeutralMax = 0.0;
        scalar nParticleIonNeutalMax = 0.0;
        forAll(neutralSpecies,i)//Species 1
        {
            const DynamicList<typename CloudType::parcelType*>& cellParcels(sortedCellOccupancy[celli][neutralSpecies[i]]);
            forAll(cellParcels,i) {
                typename CloudType::parcelType* p = cellParcels[i];
                collisionList.append(p);

                scalar nParticle = p->nParticle();
                if(nParticle > nParticleNeutralMax)
                    nParticleNeutralMax = nParticle;
            }

        }
        nParticleIonNeutalMax = nParticleNeutralMax;//If there is no higher weight for the ions use the one for the neutrals
        forAll(chargedSpecies,i)//Species 1
        {
            if(chargedSpecies[i] == cloud.electronTypeId())
                continue;

            const DynamicList<typename CloudType::parcelType*>& cellParcels(sortedCellOccupancy[celli][chargedSpecies[i]]);
            forAll(cellParcels,i) {
                typename CloudType::parcelType* p = cellParcels[i];
                ionCollisionList.append(p);

                scalar nParticle = p->nParticle();
                if(nParticle > nParticleIonNeutalMax)
                    nParticleIonNeutalMax = nParticle;
            }
        }

        //Number of particles
        nC = collisionList.size();
        iC = ionCollisionList.size();

        // -------------------------------------------------------------------------------//
        //                              Neutral - Neutral                                 //
        // -------------------------------------------------------------------------------//

        if(nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(nC);
            const point& cC = mesh.cellCentres()[celli];

            forAll(collisionList, i)
            {
                const typename CloudType::parcelType& p = *collisionList[i];
                vector relPos = p.position() - cC;

                label subCell =
                    pos0(relPos.x()) + 2*pos0(relPos.y()) + 4*pos0(relPos.z());

                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            //Calculate number of pairs that can collide
            scalar selectedPairs =
                    collisionSelectionRemainder_[celli]
                    + 0.5*nC*(nC - 1)*nParticleNeutralMax*sigmaTcRMax*deltaT
                    /mesh.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;
            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen.sampleAB<label>(0, nC);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];
                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen.sampleAB<label>(0, nSC)];
                    } while (candidateP == candidateQ);
                    }
                    else
                    {
                        // Select a possible second collision candidate from the
                        // whole cell.  If the same candidate is chosen, choose
                        // again.

                        do
                        {
                            candidateQ = rndGen.sampleAB<label>(0, nC);
                        } while (candidateP == candidateQ);
                        }

                     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     // uniform candidate selection procedure

                     // // Select the first collision candidate
                     // label candidateP = rndGen_.sampleAB<label>(0, nC);

                     // // Select a possible second collision candidate
                     // label candidateQ = rndGen_.sampleAB<label>(0, nC);

                     // // If the same candidate is chosen, choose again
                     // while (candidateP == candidateQ)
                     // {
                     //     candidateQ = rndGen_.sampleAB<label>(0, nC);
                     // }

                     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                     typename CloudType::parcelType* parcelP = collisionList[candidateP];
                     typename CloudType::parcelType* parcelQ = collisionList[candidateQ];

                     scalar sTcR = totalCrossSection_->sigmaTcR
                     (
                          *parcelP,
                          *parcelQ
                     );

                     // Update the maximum value of sigmaTcR stored, but use the
                     // initial value in the acceptance-rejection criteria because
                     // the number of collision candidates selected was based on this

                     if (sTcR > sigmaTcRMax_[celli])
                     {
                         sigmaTcRMax_[celli] = sTcR;
                     }

                     if ((sTcR/sigmaTcRMax) > rndGen.scalar01())
                     {
                         elasticCollision
                         (
                             parcelP,
                             parcelQ
                         );

                         collisions++;
                    }
             }

        }
        // -------------------------------------------------------------------------------//
        //                                  Ion - Neutral                                 //
        // -------------------------------------------------------------------------------//
        if(nC > 0 && iC > 0)
        {
            scalar sigmaTcRMax = sigmaTIoncRMax_[celli];

            //Calculate number of possible pairs
            scalar selectedPairs =
                    collisionSelectionIonRemainder_[celli]
                    + nC*iC*nParticleIonNeutalMax*sigmaTcRMax*deltaT
                    /mesh.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionSelectionIonRemainder_[celli] = selectedPairs - nCandidates;
            ionCollisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {

                //Pick random canidates
                label candidateP = rndGen.sampleAB<label>(0, nC);
                label candidateQ = rndGen.sampleAB<label>(0, iC);

                typename CloudType::parcelType* parcelP = collisionList[candidateP];
                typename CloudType::parcelType* parcelQ = ionCollisionList[candidateQ];

                scalar sTcR = totalCrossSection_->sigmaTcR
                (
                     *parcelP,
                     *parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sTcR > sigmaTIoncRMax_[celli])
                {
                    sigmaTIoncRMax_[celli] = sTcR;
                }
                scalar QeCR = sTcR;
                scalar QexcCR = 0.0;

                //Nanbu 0.5 sigmaT is the charge exchange cross section: https://doi.org/10.1109/27.887765
                if(chargeExchangeCollision_) {
                    QeCR *= 0.5;
                    QexcCR = QeCR;
                }

                scalar rndU = rndGen.scalar01();
                if((QeCR/sigmaTcRMax) > rndU)//Elastic collision
                {
                    elasticCollision
                    (
                        parcelP,
                        parcelQ
                    );
                    ionCollisions++;
                }
                else if((QeCR+QexcCR)/sigmaTcRMax > rndU)//Charge exchange collision
                {
                    chargeExchangeCollision
                    (
                        parcelP,
                        parcelQ
                    );
                    chargeExchanges++;
                    ionCollisions++;
                }
            }
        }
    }
    //Parallel COM number of events
    reduce(chargeExchanges,sumOp<label>());
    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    reduce(ionCollisions, sumOp<label>());

    reduce(ionCollisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();
    sigmaTIoncRMax_.correctBoundaryConditions();

    //Print some info
    if (collisionCandidates)
    {
        Info<< "    Collisions between neutral species  = "
                << collisions << nl
                << "    Acceptance rate                     = "
                << scalar(collisions)/scalar(collisionCandidates) << nl
                << nl << endl;
    }
    else
    {
        Info<< "    No collisions between neutral species" << nl << endl;
    }

    if (ionCollisionCandidates)
    {
        Info<< "    Collisions between ion/neutral species = "
                << ionCollisions << nl
                << "    Elastic                                = "
                << ionCollisions-chargeExchanges << nl
                << "    Charge exchange                        = "
                << chargeExchanges << nl
                << "    Acceptance rate                        = "
                << scalar(ionCollisions)/scalar(ionCollisionCandidates) << nl
                << nl << endl;
    }
    else
    {
        Info<< "    No collisions between ion/neutral species" << nl << endl;
    }
}
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::performBackgroundCollisions()
{

    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    const List<label>& ionSpecies(cloud.ionSpecies());
    scalar dt = mesh.time().deltaTValue();


    DynamicList<typename CloudType::parcelType*> collisionList;
    DynamicList<typename CloudType::parcelType*> ionCollisionList;

    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    //Used for counting events
    label collisionsIon = 0;
    label collisionsNeutral = 0;
    label collisionCandidatesIon = 0;
    label collisionCandidatesNeutral = 0;
    label chargeExchanges = 0;

    //update background (can have changed) FIXME: use update flag?
    tmp<scalarField> tbgNumberDensity = backgroundGas.numberDensity();
    backgroundCollisionPropIon_.field() = 1.0 -exp(-sigmaTIoncRMax_*tbgNumberDensity*dt);
    backgroundCollisionProp_.field() = 1.0 -exp(-sigmaTcRMax_*tbgNumberDensity*dt);

    forAll(sortedCellOccupancy, celli)// Go through all cells and add the particles to the collision lists
    {
        label nN = 0;
        label nI = 0;
        collisionList.clear();
        ionCollisionList.clear();

        forAll(neutralSpecies,i)
        {
            label typeId = neutralSpecies[i];
            const DynamicList<typename CloudType::parcelType*>& cellParcels(sortedCellOccupancy[celli][typeId]);

            forAll(cellParcels,j)
            {
                typename CloudType::parcelType* p = cellParcels[j];
                collisionList.append(p);
                nN++;
            }
        }
        forAll(ionSpecies,i)
        {
            label typeId = ionSpecies[i];
            const DynamicList<typename CloudType::parcelType*>& cellParcels(sortedCellOccupancy[celli][typeId]);
            forAll(cellParcels,j)
            {
                typename CloudType::parcelType* p = cellParcels[j];
                ionCollisionList.append(p);
                nI++;
            }
        }

        //Collision of neutrals with the background
        if(nN > 0)
        {
            scalar collisions = backgroundCollisionProp_[celli]*nN+backgroundCSR_[celli];//Calculate number of collisions
            label nCandidates(collisions);
            backgroundCSR_[celli] = collisions-nCandidates;
            collisionCandidatesNeutral += nCandidates;
            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            for (label i = 0; i < nCandidates; i++)
            {
                //Pick random collision partners
                label candidate = rndGen.sampleAB<label>(0, nN);//FIXME: candidate should only be selected ONCE !!!!!
                typename CloudType::parcelType* parcel = collisionList[candidate];

                vector U = backgroundGas.sampleVelocity(celli);
                vector preU = U;
                label bgTypeId = backgroundGas.species();
                scalar sTcR = totalCrossSection_->sigmaTcR
                (
                     *parcel,
                     U,
                     bgTypeId
                );

                //Update the maximum value
                if (sTcR > sigmaTcRMax_[celli])
                {
                    sigmaTcRMax_[celli] = sTcR;
                    backgroundCollisionProp_[celli] = 1.0 -::exp(-sTcR*backgroundGas.numberDensity().operator ()()[celli]*dt);
                }

                if ((sTcR/sigmaTcRMax) > rndGen.scalar01())
                {
                    collisionsNeutral++;
                    elasticCollision//Perform the collision
                    (
                        parcel,
                        U,
                        bgTypeId
                    );
                    backgroundGas.collisionUpdate(CollisionEvent::Elastic, celli, preU, U);
                }
            }

        }
        if(nI > 0)//Collisions of ions with the background
        {
            scalar collisions = backgroundCollisionPropIon_[celli]*nI+backgroundCSRIon_[celli];//Number of collisions
            label nCandidates(collisions);
            collisionCandidatesIon += nCandidates;
            backgroundCSRIon_[celli] = collisions-nCandidates;

            scalar sigmaTcRMax = sigmaTIoncRMax_[celli];

            for (label i = 0; i < nCandidates; i++)
            {
                //Pick random collision partners
                label candidate = rndGen.sampleAB<label>(0, nI);//FIXME: candidate should only be selected ONCE !!!!!
                typename CloudType::parcelType* parcel = ionCollisionList[candidate];

                vector U = backgroundGas.sampleVelocity(celli);
                vector preU = U;
                label bgTypeId = backgroundGas.species();
                scalar sTcR = totalCrossSection_->sigmaTcR
                (
                     *parcel,
                     U,
                     bgTypeId
                );

                //Update the maximum value
                if (sTcR > sigmaTIoncRMax_[celli])
                {
                    sigmaTIoncRMax_[celli] = sTcR;
                    backgroundCollisionPropIon_[celli] = 1.0 -::exp(-sTcR*backgroundGas.numberDensity().operator ()()[celli]*dt);
                }

                scalar QeCR = sTcR;
                scalar QexcCR = 0.0;

                //Nanbu 0.5 sigmaT is the charge exchange cross section: https://doi.org/10.1109/27.887765
                if(chargeExchangeCollision_) {
                    QeCR *= 0.5;
                    QexcCR = QeCR;
                }

                scalar rndU = rndGen.scalar01();
                if((QeCR/sigmaTcRMax) > rndU)//Elastic collision
                {
                    elasticCollision
                    (
                        parcel,
                        U,
                        bgTypeId
                    );
                    backgroundGas.collisionUpdate(CollisionEvent::Elastic, celli, preU, U);
                    collisionsIon++;
                }
                else if((QeCR+QexcCR)/sigmaTcRMax > rndU)//Charge exchange
                {
                    chargeExchangeCollision
                    (
                        parcel,
                        U,
                        bgTypeId
                    );
                    backgroundGas.collisionUpdate(CollisionEvent::ChargeExchange, celli, preU, U);

                    chargeExchanges++;
                    collisionsIon++;
                }
            }
        }

    }
    sigmaTcRMax_.correctBoundaryConditions();
    sigmaTIoncRMax_.correctBoundaryConditions();

    //Parallel COM number of events
    reduce(chargeExchanges,sumOp<label>());
    reduce(collisionsIon, sumOp<label>());
    reduce(collisionsNeutral, sumOp<label>());

    reduce(collisionCandidatesIon, sumOp<label>());
    reduce(collisionCandidatesNeutral, sumOp<label>());

    //Print info
    if (collisionCandidatesNeutral)
    {
        Info << "    Collisions between neutals and background  = "
             << collisionsNeutral << nl
             << "    Acceptance rate                                      = "
             << scalar(collisionsNeutral)/scalar(collisionCandidatesNeutral) << nl
             << nl << endl;
    }
    else
        Info<< "    No Collisions between neutals and the background" << nl << endl;


    if (collisionCandidatesIon)
    {
        Info << "    Collisions between ions and background  = "
             << collisionsIon << nl
             << "    Elastic                                = "
             << collisionsIon-chargeExchanges << nl
             << "    Charge exchange                        = "
             << chargeExchanges << nl
             << "    Acceptance rate                                      = "
             << scalar(collisionsIon)/scalar(collisionCandidatesIon) << nl
             << nl << endl;
    }
    else
        Info<< "    No Collisions between ions and the background" << nl << endl;

}

template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::handleCollisions()
{
    Info << "[Binary Collisions]" << nl;
    performCollisions();

    const BackgroundGasModel<CloudType>& backgroundGas(this->owner().backgroundGas());
    if(backgroundGas.active())
        performBackgroundCollisions();
}

/*
Foam::BinaryCollisionModel<CloudType>::initialize called by picInitialise
*/
template<class CloudType>
void Foam::BinaryCollisionModel<CloudType>::initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities)
{
    CloudType& cloud(this->owner());
    const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    if(neutralSpecies.size() == 0)
    {
        WarningInFunction << "No neutral species defined!" << endl;
        return;
    }


    label maxSpecies = neutralSpecies[0];//use first species if none is set
    scalar maxDensity = 0.0;

    //Look for the most common neutral species
    scalar maxbgDensity = max(backgroundGas.numberDensity());
    scalar maxbgTemperature = max(backgroundGas.temperature());
    forAll(neutralSpecies,i)
    {
        label typeId = neutralSpecies[i];
        scalar curDensity = numberDensities[typeId];

        if(backgroundGas.active() && backgroundGas.species() == typeId)
            curDensity += maxbgDensity;

        if(curDensity > maxDensity)
        {
            maxDensity = curDensity;
            maxSpecies = typeId;
        }
    }

    const typename CloudType::parcelType::constantProperties& cP(cloud.constProps(maxSpecies));

    scalar T = temperatures[maxSpecies];

    if(backgroundGas.active() && maxSpecies == backgroundGas.species())
    {
        if(maxbgTemperature > T)
            T = maxbgTemperature;
    }

    if(T < VSMALL) {
        Warning << "    In initialization of sigmaTcRMax_: Zero temperature for " << cloud.typeIdList()[maxSpecies] << " species" << nl << "    |_ Using default value of 293.15 K..." << endl;
        T = 293.15;//Default value was chosen to be 293.15
    }

    //Initialize the sigmaTcRMax field
    sigmaTcRMax_.primitiveFieldRef() = cP.sigmaT()*cloud.maxwellianMostProbableSpeed
    (
        T,
        cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();

    const List<label>& ionSpecies(cloud.ionSpecies());
    if(ionSpecies.size() == 0)
        return;

    //Look for the most common ion species
    label maxIonSpecies = -1;
    scalar maxIonDensity = 0.0;
    forAll(ionSpecies,i)
    {
        label typeId = ionSpecies[i];

        if(numberDensities[typeId] > maxIonDensity)
        {
            maxIonDensity = numberDensities[typeId];
            maxIonSpecies = typeId;
        }
    }

    if(maxIonSpecies == -1)//no inital ions use neutral max
    {
        sigmaTIoncRMax_.primitiveFieldRef() = cP.sigmaT()*cloud.maxwellianMostProbableSpeed
        (
              T,
              cP.mass()
        );
    }
    else
    {
        const typename CloudType::parcelType::constantProperties& cPIon = cloud.constProps(maxIonSpecies);


        scalar Tion = temperatures[maxIonSpecies];
        if(Tion < VSMALL) {
            Warning << "    In initialization of sigmaTIoncRMax_: Zero temperature for ion species" << nl << "    |_ Using default value of 1 eV..." << endl;
            Tion = 11604.52;//Default value was chosen to be 1 eV
        }

        scalar Ui = cloud.maxwellianMostProbableSpeed
        (
             Tion,
             cPIon.mass()
        );


        //Use the hard sphere model and update the sigmaTIoncRMax field
        scalar dPQ = 0.5*(cPIon.d()+cP.d());
        sigmaTIoncRMax_.primitiveFieldRef() = dPQ*dPQ*constant::mathematical::pi*Ui;
    }

    sigmaTIoncRMax_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BinaryCollisionModelNew.C"

// ************************************************************************* //
