/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 picFoam
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

#include "IonNeutralCollisionModel.H"
#include "BackgroundGasModel.H"
#include "Random.H"
#include "constants.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IonNeutralCollisionModel<CloudType>::IonNeutralCollisionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    elasticCSModels_(owner),
    chargeExCSModels_(owner),
    backgroundSigmaTcRMax_(
        IOobject
        (
            "picSigmaTcRMax_IonBackground",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.mesh(),
        dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), 0.0)
    ),
    backgroundCollisionRemainder_(
        IOobject
        (
            "backgroundCollisionRemainder_IonNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_IonNeutral", dimless, 0)
   ),
   sigmaTcRMax_
   (
      IOobject
      (
          "picSigmaTcRMax_IonNeutral",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), 0.0)
  ),
  collisionRemainder_
  (
      IOobject
      (
          "collisionRemainder_IonNeutral",
          owner.mesh().time().timeName(),
          owner.mesh()
      ),
      owner.mesh(),
      dimensionedScalar("collisionRemainder_IonNeutral", dimless, 0)
  ),
  weightCorrection_(),
  updateNeutralVelocity_(true)
{}


template<class CloudType>
Foam::IonNeutralCollisionModel<CloudType>::IonNeutralCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict("IonNeutralCollisionCoeffs")),
    elasticCSModels_
    (
       coeffDict_.subDict("CrossSectionModels"),
       owner
    ),
    chargeExCSModels_
    (
       coeffDict_.subDict("CrossSectionModels"),
       owner
    ),
    backgroundSigmaTcRMax_(
        IOobject
        (
            "picSigmaTcRMax_IonBackground",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        owner.mesh(),
        dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    backgroundCollisionRemainder_(
        IOobject
        (
            "backgroundCollisionRemainder_IonNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_IonNeutral", dimless, 0)
   ),
   sigmaTcRMax_
   (
      IOobject
      (
          "picSigmaTcRMax_IonNeutral",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), Zero),
      calculatedFvPatchScalarField::typeName
  ),
  collisionRemainder_
  (
      IOobject
      (
          "collisionRemainder_IonNeutral",
          owner.mesh().time().timeName(),
          owner.mesh()
      ),
      owner.mesh(),
      dimensionedScalar("collisionRemainder_IonNeutral", dimless, 0)
  ),
  weightCorrection_
  (
      WeightCorrectionModel<CloudType>::New
      (
          coeffDict_,
          owner
      )
  ),
  updateNeutralVelocity_(coeffDict_.lookupOrDefault("updateNeutralVelocity",true))
{
    if(owner.neutralSpecies().size() > 1 || owner.ionSpecies().size() > 1)
    {
        //FIXME: Model only supports one species at the moment, because we reuse the CrossSection model of the ElectronNeutral model
        FatalErrorInFunction << "This model is designed to support one neutral and one ion species only!" << abort(FatalError);
    }

    //Read above, will be zero if no file was found...
    Info << "|->    Reading picSigmaTcRMax_IonNeutral if present..." << nl
         << "           sigmaTcRmax: " << gMax(sigmaTcRMax_.primitiveFieldRef()) << endl;

    //Print info on the background gas model
    const BackgroundGasModel<CloudType>& backgroundGas(owner.backgroundGas());
    if(backgroundGas.active())
    {
        label bgTypeId = backgroundGas.species();

        label ionTypeId = owner.ionTypeId(bgTypeId);
        if(ionTypeId < 0)
            FatalErrorInFunction << "No ion type for neutral species " << owner.typeIdList()[bgTypeId] << abort(FatalError);

        Info << "|->    BackgroundGasModel " << backgroundGas.type() << " ion collision properties:" << nl
             << "           species: " << owner.typeIdList()[bgTypeId] << nl
             << "           avg number density: " << gAverage(backgroundGas.numberDensity()) << " m^-3" << nl
             << "           avg temperature: " << gAverage(backgroundGas.temperature()) << " K" << nl
             << "           sigmaTcRmax: " << gMax(backgroundSigmaTcRMax_) << nl << endl;
    }

    forAll(collisionRemainder_, i)
    {
        backgroundCollisionRemainder_[i] = owner.rndGen().scalar01();
        collisionRemainder_[i] = owner.rndGen().scalar01();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IonNeutralCollisionModel<CloudType>::~IonNeutralCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::IonNeutralCollisionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType&
Foam::IonNeutralCollisionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary&
Foam::IonNeutralCollisionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::IonNeutralCollisionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::handleCollisions()
{
    Info << "[Ion-Neutral Collisions]" << nl;
    performCollisions();

    const BackgroundGasModel<CloudType>& backgroundGas(this->owner().backgroundGas());
    if(backgroundGas.active())
        performBackgroundCollisions();
}

template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities)
{
    if(!active())
        return;

    Info << "IonNeutral CollisionModel:" << endl;

    CloudType& cloud(this->owner());
    const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    if(neutralSpecies.size() == 0)
    {
        Info << "|->    WARNING IonNeutralCollisionModel::initialize() no neutral species defined!" << endl;
        return;
    }

    //Lookup initialization mode
    dictionary& crossSectionDict(cloud.picInitialiseDict().subDict("CrossSectionInitialization"));
    word initMode = crossSectionDict.lookup("ionNeutral");
    if(initMode != "FindMaximum" && initMode != "Temperature")
        FatalErrorInFunction << "Invalid initialization mode, valid options are: (FindMaximum, Temperature)" << abort(FatalError);

    bool findMaxium = false;
    if(initMode == "FindMaximum")
        findMaxium = true;

    //FindMaxumum mode: search for the maximum cross section in given range
    if(findMaxium)
    {
        //Read paramater from picInitialiseDict
        dictionary& searchParamDict(crossSectionDict.subDict("ionNeutral_searchParameter"));
        scalar maximumEnergy = readScalar(searchParamDict.lookup("maximumEnergy"));
        scalar energyStep = readScalar(searchParamDict.lookup("energyStep"));

        Info << "|->    Looking up and initializing sigmaTcRMax to the maximum value in range 0.0 to " << maximumEnergy << " eV..." << endl;
        forAll(neutralSpecies,i)//Note list has only one entry...
        {
            label typeId = neutralSpecies[i];
            label ionTypeId = cloud.ionSpecies()[typeId];
            const typename CloudType::parcelType::constantProperties& cP(cloud.constProps(ionTypeId));
            scalar massI = cP.mass();

            //Search maximum
            scalar eV = 0.0;
            scalar sigmaTMax = 0.0;
            for(;eV < maximumEnergy; eV+=energyStep)
            {
                scalar Qel = elasticCrossSections()[typeId].crossSection(eV);
                scalar Qex = chargeExCrossSections()[typeId].crossSection(eV);

                sigmaTMax = max(sigmaTMax,::sqrt(2.0*constant::electromagnetic::e.value()*eV/massI)*(Qel+Qex));
            }
            //Store maximum
            sigmaTcRMax_.primitiveFieldRef() = sigmaTMax;
        }

        //Background gas model initialization
        const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());
        if(backgroundGas.active())
        {
            label bgSpecies = backgroundGas.species();
            label bgIonSpecies = cloud.ionTypeId(bgSpecies);//-1 checked in constructor

            const typename CloudType::parcelType::constantProperties& cPIon(cloud.constProps(bgIonSpecies));
            scalar massI = cPIon.mass();

            scalar eV = 0.0;
            scalar sigmaTMax = 0.0;
            for(;eV < maximumEnergy; eV+=energyStep)
            {
                scalar Qel = elasticCrossSections()[bgSpecies].crossSection(eV);
                scalar Qex = chargeExCrossSections()[bgSpecies].crossSection(eV);

                sigmaTMax = max(sigmaTMax,::sqrt(2.0*constant::electromagnetic::e.value()*eV/massI)*(Qel+Qex));
            }
            backgroundSigmaTcRMax_.primitiveFieldRef() = sigmaTMax;
        }
    }
    else//Temperature mode: Use specified temperature for initialization
    {
        Info << "|->    Initializing sigmaTcRMax by specified temperatures..." << endl;
        forAll(neutralSpecies,i)//Note list has only one entry...
        {
            label typeId = neutralSpecies[i];
            label ionTypeId = cloud.ionSpecies()[typeId];

            scalar Tion = temperatures[ionTypeId];
            if(Tion < VSMALL) {
                Info << "|->    WARNING in initialization of picSigmaTcRMax_IonNeutral: Zero temperature for ion species" << nl << "       |= Using default value of 1 eV..." << endl;
                Tion = 11604.52;
            }

            const typename CloudType::parcelType::constantProperties& cP(cloud.constProps(ionTypeId));
            scalar massI = cP.mass();

            scalar eV = Tion*constant::physicoChemical::k.value()/constant::electromagnetic::e.value();

            scalar Qel = elasticCrossSections()[typeId].crossSection(eV);
            scalar Qex = chargeExCrossSections()[typeId].crossSection(eV);
            scalar Qt = Qel+Qex;

            sigmaTcRMax_.primitiveFieldRef() = Qt*cloud.maxwellianMostProbableSpeed
            (
                 Tion,
                 massI
            );

            sigmaTcRMax_.correctBoundaryConditions();

        }

        //Background gas model initialization
        if(backgroundGas.active())
        {
            label bgSpecies = backgroundGas.species();

            //Right now: bgIonSpecies == ionTypeId and Tion == T, this is here for future extensions... performance doesn't matter during initialization
            label bgIonSpecies = cloud.ionTypeId(bgSpecies);//-1 checked in constructor
            scalar Tion = temperatures[bgIonSpecies];
            if(Tion < VSMALL) {
                Info << "|->    WARNING in initialization of picSigmaTcRMax_IonBackground: Zero temperature for ion species" << nl << "       |= Using default value of 1 eV..." << endl;
                Tion = 11604.52;
            }

            const typename CloudType::parcelType::constantProperties& cPIon(cloud.constProps(bgIonSpecies));
            scalar eV = Tion*constant::physicoChemical::k.value()/constant::electromagnetic::e.value();

            scalar Qel = elasticCrossSections()[bgSpecies].crossSection(eV);
            scalar Qex = chargeExCrossSections()[bgSpecies].crossSection(eV);
            scalar Qt = Qel+Qex;

            backgroundSigmaTcRMax_.primitiveFieldRef() = Qt*cloud.maxwellianMostProbableSpeed
            (
                 Tion,
                 cPIon.mass()
            );
            backgroundSigmaTcRMax_.correctBoundaryConditions();
        }
    }
}

//Charge exchange collision particle-particle
template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::chargeExchangeCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ)
{
    scalar rndU = this->owner().rndGen().scalar01();

    //Identity Switch
    vector tmpP = pP->U();
    if(rndU <= pQ->nParticle()/pP->nParticle())
        pP->U() = pQ->U();

    if(rndU <= pP->nParticle()/pQ->nParticle())
        pQ->U() = tmpP;

    if(!this->updateNeutralVelocity())
        pP->U() = tmpP;
}

//Charge exchange collision background gas
template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::chargeExchangeCollision(typename CloudType::parcelType* pP, vector& Uq)
{
    vector tmpP = pP->U();
    pP->U() = Uq;
    Uq = tmpP;
}

//Elastic collision background gas
template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::elasticCollision(typename CloudType::parcelType* pP, vector& Uq)
{
    BackgroundGasModel<CloudType>& backgroundGas(this->owner().backgroundGas());
    label typeId = backgroundGas.species();

    updateVelocity(*pP,Uq,typeId);
}

//Elastic collision particle-particle
template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::elasticCollision(typename CloudType::parcelType* pP, typename CloudType::parcelType* pQ)
{
    vector preUp = pP->U();
    vector preUq = pQ->U();

    updateVelocity(*pP,*pQ);
    weightCorrection_->correctVelocity(pP,pQ,preUp,preUq);


    if(!this->updateNeutralVelocity())
        pP->U() = preUp;
}

//Handle background gas collisions
template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::performBackgroundCollisions()
{

    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());
    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label bgTypeId = backgroundGas.species();
    label ionTypeId = cloud.ionTypeId(bgTypeId);

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

    label collisionsElastic = 0;
    label collisionsChargeEx = 0;
    label collisionCandidatesIon = 0;

    //Used for most probable speed...
    //const scalarField& bgTemperature(backgroundGas.temperature());
    //scalar massN = cloud.constProps(backgroundGas.species()).mass();

    scalar dt = mesh.time().deltaTValue();
    scalar massI = cloud.constProps(ionTypeId).mass();
    const scalarField& bgDensity(backgroundGas.numberDensity());

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        //Number of ions
        const DynamicList<typename CloudType::parcelType*>& ionParcels(sortedCellOccupancy[celli][ionTypeId]);
        label nI = ionParcels.size();

        scalar sigmaTcRMax = backgroundSigmaTcRMax_[celli];

        //Calculate the number of candidates for a collision
        scalar collisions = backgroundCollisionRemainder_[celli] + nI*bgDensity[celli]*sigmaTcRMax*dt;

        //Proposed in Nanbu2000, used in xpdp1, not in oopd1...
        //scalar collisions = backgroundCollisionRemainder_[celli] + (1.0 - exp(-sigmaTcRMax*bgDensity[celli]*dt))*nI;

        label nCandidates(collisions);
        collisionCandidatesIon += nCandidates;

        //Store remainder
        backgroundCollisionRemainder_[celli] = collisions-nCandidates;

        //Go through nCandidates and check if a collision occurs
        for (label i = 0; i < nCandidates; i++)
        {
            //Pick a random ion
            label candidate = rndGen.sampleAB<label>(0, nI);//FIXME: candidate should only be selected ONCE !!!!!
            typename CloudType::parcelType* parcel = ionParcels[candidate];

            //Sample velocity / Alternative: most probable speed = sqrt(2.0*constant::physicoChemical::k.value()*bgTemperature[celli]/massN)
            vector Un = backgroundGas.sampleVelocity(celli);
            vector preUn = Un;

            scalar cR = mag(parcel->U()-Un);
            scalar gamma=1.0/::sqrt(1.0-((cR*cR)/(constant::universal::c.value()*constant::universal::c.value())));

            //Calc relativistic kin energy
            scalar eVkinEnergy = (gamma-1.0)*massI*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());

            scalar sigmaT = sigmaTcRMax/cR;

            //Calculate the cross section
            scalar Qel = elasticCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qex = chargeExCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qt = Qel+Qex;

            scalar exThreshold = chargeExCrossSections()[bgTypeId].threshold();

            //Update backgroundSigmaTcRMax_... this needs to be done because of fixed/constant cross sections even when FindMaximum mode is used cR can make this larger...
            if(Qt*cR > backgroundSigmaTcRMax_[celli]) {
                backgroundSigmaTcRMax_[celli] = Qt*cR;
            }

            //Compare against the maximum cross section
            scalar rndU = rndGen.scalar01();
            if(rndU <= Qel/sigmaT)
            {
                elasticCollision
                (
                    parcel,
                    Un
                );
                backgroundGas.collisionUpdate(CollisionEvent::Elastic, celli, preUn, Un);
                collisionsElastic++;
            }
            else if(eVkinEnergy >= exThreshold && rndU <= Qt/sigmaT)
            {
                chargeExchangeCollision
                (
                    parcel,
                    Un
                );
                backgroundGas.collisionUpdate(CollisionEvent::ChargeExchange, celli, preUn, Un);
                collisionsChargeEx++;
            }
        }
    }


    //Print some information
    reduce(collisionsElastic,sumOp<label>());
    reduce(collisionsChargeEx, sumOp<label>());
    reduce(collisionCandidatesIon, sumOp<label>());

    if (collisionCandidatesIon)
    {
        Info << "    Collisions between ions and background  = "
             << collisionsElastic+collisionsChargeEx << nl
             << "    Elastic                                = "
             << collisionsElastic << nl
             << "    Charge exchange                        = "
             << collisionsChargeEx << nl
             << "    Acceptance rate                                      = "
             << scalar(collisionsElastic+collisionsChargeEx)/scalar(collisionCandidatesIon) << nl
             << nl << endl;
    }
    else
        Info<< "    No Collisions between ions and the background" << nl << endl;
}


template<class CloudType>
void Foam::IonNeutralCollisionModel<CloudType>::performCollisions()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());

    const List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    const List<label>& chargedSpecies(cloud.chargedSpecies());

    label ionCollisionCandidates = 0;
    label ionCollisions = 0;
    label chargeExchanges = 0;

    DynamicList<typename CloudType::parcelType*> collisionList;
    DynamicList<typename CloudType::parcelType*> ionCollisionList;

    scalar dt = mesh.time().deltaTValue();

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        label nC(0), iC(0);
        collisionList.clear();
        ionCollisionList.clear();
        scalar nParticleNeutralMax = 0.0;
        scalar nParticleIonNeutalMax = 0.0;

        //Get the maximum neutral parcel weight in this cell
        forAll(neutralSpecies,i)//Note contains only one species
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

        //Get the maximum ion parcel weight in this cell
        nParticleIonNeutalMax = nParticleNeutralMax;
        forAll(chargedSpecies,i)//Note contains only one species
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

        nC = collisionList.size();
        iC = ionCollisionList.size();

        // -------------------------------------------------------------------------------//
        //                                  Ion - Neutral                                 //
        // -------------------------------------------------------------------------------//
        if(nC > 0 && iC > 0)
        {
            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            //Calculate number of candidates
            scalar selectedPairs =
                    collisionRemainder_[celli]
                    + nC*iC*nParticleIonNeutalMax*sigmaTcRMax*dt
                    /mesh.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionRemainder_[celli] = selectedPairs - nCandidates;
            ionCollisionCandidates += nCandidates;

            //Check nCandidates
            for (label c = 0; c < nCandidates; c++)
            {
                //Pick random parcels
                label candidateP = rndGen.sampleAB<label>(0, nC);
                label candidateQ = rndGen.sampleAB<label>(0, iC);

                typename CloudType::parcelType* parcelP = collisionList[candidateP];
                label neutralTypeId = parcelP->typeId();
                typename CloudType::parcelType* parcelQ = ionCollisionList[candidateQ];
                label ionTypeId = parcelQ->typeId();


                scalar massI = cloud.constProps(ionTypeId).mass();
                //Relative velocity
                scalar cR = mag(parcelQ->U()-parcelP->U());
                scalar gamma=1.0/::sqrt(1.0-((cR*cR)/(constant::universal::c.value()*constant::universal::c.value())));

                //Calc relativistic kin energy
                scalar eVkinEnergy = (gamma-1.0)*massI*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());

                //Calculate total cross section
                scalar Qel = elasticCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                scalar Qex = chargeExCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                scalar Qt = Qel+Qex;

                scalar exThreshold = chargeExCrossSections()[neutralTypeId].threshold();

                //Update sigmaTcRMax_ if the current one is large
                if(Qt*cR > sigmaTcRMax_[celli]) {
                    sigmaTcRMax_[celli] = Qt*mag(cR);
                }

                //Check if a collision occurs
                scalar rndU = rndGen.scalar01();
                if(rndU <= Qel*cR/sigmaTcRMax)
                {
                    ionCollisions++;
                    elasticCollision(parcelP, parcelQ);
                }
                else if(eVkinEnergy >= exThreshold && rndU <= Qt*cR/sigmaTcRMax)
                {
                    chargeExchanges++;
                    chargeExchangeCollision(parcelP,parcelQ);
                }
            }
        }
    }

    //Print some information
    reduce(chargeExchanges,sumOp<label>());
    reduce(ionCollisions, sumOp<label>());
    reduce(ionCollisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (ionCollisionCandidates)
    {
        Info<< "    Collisions between ion/neutral species = "
                << ionCollisions+chargeExchanges << nl
                << "    Elastic                                = "
                << ionCollisions << nl
                << "    Charge exchange                        = "
                << chargeExchanges << nl
                << "    Acceptance rate                        = "
                << scalar(ionCollisions+chargeExchanges)/scalar(ionCollisionCandidates) << nl
                << nl << endl;
    }
    else
    {
        Info<< "    No collisions between ion/neutral species" << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IonNeutralCollisionModelNew.C"

// ************************************************************************* //
