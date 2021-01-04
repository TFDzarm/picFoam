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

#include "ElectronNeutralCollisionModel.H"
#include "BackgroundGasModel.H"
#include "Random.H"
#include "constants.H"
#include "NanbuCorrection.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectronNeutralCollisionModel<CloudType>::ElectronNeutralCollisionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    elasticCSModels_(owner),
    excitationCSModels_(owner),
    ionizationCSModels_(owner),
    deleteNeutral_(true),
    updateNeutralVelocity_(true),
    ionizationNEquivalentParticles_(0.0),
    backgroundSigmaTcRMax_(0.0),
    backgroundCollisionProp_(owner.mesh().cells().size(),0.0),
    backgroundCollisionRemainder_(
        IOobject
        (
            "backgroundCollisionRemainder_ElectronNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_ElectronNeutral", dimless, 0)
   ),
   sigmaTcRMax_
   (
      IOobject
      (
          "picSigmaTcRMax_ElectronNeutral",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      owner.mesh()
  ),
  collisionRemainder_
  (
      IOobject
      (
          "collisionRemainder_ElectronNeutral",
          owner.mesh().time().timeName(),
          owner.mesh()
      ),
      owner.mesh(),
      dimensionedScalar("collisionRemainder_ElectronNeutral", dimless, 0)
  ),
  weightCorrection_(),
  ionizationAverage_(0.0),
  averageStart_(owner.mesh().time().value()),
  elasticCollisions_
  (
      IOobject
      (
          "elasticElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero", dimless, 0.0),
      calculatedFvPatchScalarField::typeName
  ),
  excitationCollisions_
  (
      IOobject
      (
          "excitationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero", dimless, 0.0),
      calculatedFvPatchScalarField::typeName
  ),
  ionizationCollisions_
  (
      IOobject
      (
          "ionizationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero", dimless, 0.0),
      calculatedFvPatchScalarField::typeName
  )
{}


template<class CloudType>
Foam::ElectronNeutralCollisionModel<CloudType>::ElectronNeutralCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict("ElectronNeutralCollisionCoeffs")),
    elasticCSModels_//init the elastic cross section model
    (
       coeffDict_.subDict("CrossSectionModels"),
       owner
    ),
    excitationCSModels_//init the excitation cross section model
    (
       coeffDict_.subDict("CrossSectionModels"),
       owner
    ),
    ionizationCSModels_//init the ionization cross section model
    (
       coeffDict_.subDict("CrossSectionModels"),
       owner
    ),
    deleteNeutral_(coeffDict_.lookupOrDefault("deleteNeutral",true)),
    updateNeutralVelocity_(coeffDict_.lookupOrDefault("updateNeutralVelocity",true)),
    ionizationNEquivalentParticles_(coeffDict_.lookupOrDefault("ionization_nEquivalentParticles",0.0)),
    backgroundSigmaTcRMax_(0.0),
    backgroundCollisionProp_(owner.mesh().cells().size(),0.0),
    backgroundCollisionRemainder_(
        IOobject
        (
            "backgroundCollisionRemainder_ElectronNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_ElectronNeutral", dimless, 0)
   ),
   sigmaTcRMax_
   (
      IOobject
      (
          "picSigmaTcRMax_ElectronNeutral",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      owner.mesh()
  ),
  collisionRemainder_
  (
      IOobject
      (
          "collisionRemainder_ElectronNeutral",
          owner.mesh().time().timeName(),
          owner.mesh()
      ),
      owner.mesh(),
      dimensionedScalar("collisionRemainder_ElectronNeutral", dimless, 0)
  ),
  weightCorrection_
  (
      WeightCorrectionModel<CloudType>::New
      (
          coeffDict_,
          owner
      )
  ),
  ionizationAverage_(0.0),
  averageStart_(owner.mesh().time().value()),
  elasticCollisions_
  (
      IOobject
      (
          "elasticElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
     owner.mesh()
  ),
  excitationCollisions_
  (
      IOobject
      (
          "excitationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      owner.mesh()
  ),
  ionizationCollisions_
  (
      IOobject
      (
          "ionizationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      owner.mesh()
  )
{
    const BackgroundGasModel<CloudType>& backgroundGas(owner.backgroundGas());

    //If we have a background gas calculate the maximum total cross section and collision probability
    if(backgroundGas.active())
    {
        label bgTypeId = backgroundGas.species();

        scalar dt = owner.mesh().time().deltaTValue();
        scalar massE = owner.constProps(owner.electronTypeId()).mass();
        scalar eV = 0.0;
        for(;eV < 500.0; eV+=0.1)//FIXME: use max from CS-Models? ... But for now the max is definitly << 500.0
        {
            scalar Qel = elasticCrossSections()[bgTypeId].crossSection(eV);
            scalar Qex = excitationCrossSections()[bgTypeId].crossSection(eV);
            scalar Qiz = ionizationCrossSections()[bgTypeId].crossSection(eV);

            backgroundSigmaTcRMax_ =  max(backgroundSigmaTcRMax_,::sqrt(2.0*constant::electromagnetic::e.value()*eV/massE)*(Qel+Qex+Qiz));
            //collision probability
            backgroundCollisionProp_ = 1.0 -exp(-backgroundSigmaTcRMax_*backgroundGas.numberDensity()*dt);
        }

        Info << "|->    BackgroundGasModel " << backgroundGas.type() << " electron collision properties:" << nl
             << "           species: " << owner.typeIdList()[bgTypeId] << nl
             << "           avg number density: " << gAverage(backgroundGas.numberDensity()) << " m^-3" << nl
             << "           avg temperature: " << gAverage(backgroundGas.temperature()) << " K" << nl
             << "           sigmaTcRmax: " << backgroundSigmaTcRMax_ << nl
             << "           collision Pmax: " << gMax(backgroundCollisionProp_) << nl << endl;
    }
    forAll(collisionRemainder_, i)
    {
        backgroundCollisionRemainder_[i] = owner.rndGen().scalar01();
        collisionRemainder_[i] = owner.rndGen().scalar01();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectronNeutralCollisionModel<CloudType>::~ElectronNeutralCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::ElectronNeutralCollisionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType&
Foam::ElectronNeutralCollisionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary&
Foam::ElectronNeutralCollisionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::ElectronNeutralCollisionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::handleCollisions()
{
    Info << "[Electron-Neutral Collisions]" << nl;
    performCollisions();

    const BackgroundGasModel<CloudType>& backgroundGas(this->owner().backgroundGas());
    if(backgroundGas.active())
        performBackgroundCollisions();
}

/*
Foam::ElectronNeutralCollisionModel<CloudType>::initialize called by picInitialise
Setup the total cross section
*/
template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::initialize(Field<scalar>& temperatures, Field<scalar>& numberDensities)
{
    if(!active())
        return;

    CloudType& cloud(this->owner());
    const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    if(neutralSpecies.size() == 0)
    {
        Info << "Warning ElectronNeutralCollisionModel::initialize() no neutral species defined!" << endl;
        return;
    }

    //Calculate the maximum neutral gas density
    scalar maxbgDensity = max(backgroundGas.numberDensity());
    label maxSpecies = neutralSpecies[0];//use first species if none is set
    scalar maxDensity = 0.0;
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

    scalar T = temperatures[cloud.electronTypeId()];
    if(T < VSMALL) {
        Info << "    WARNING in initialization of sigmaTElectroncRMax: Zero temperature for electron species" << nl << "    |_ Using default value of 1 eV..." << endl;
        T = 11604.52;
    }
    const typename CloudType::parcelType::constantProperties& cP(cloud.constProps(cloud.electronTypeId()));
    scalar eVkinEnergy = T*constant::physicoChemical::k.value()/constant::electromagnetic::e.value();

    //Calculate the total cross section
    scalar Qel = elasticCrossSections()[maxSpecies].crossSection(eVkinEnergy);
    scalar Qex = excitationCrossSections()[maxSpecies].crossSection(eVkinEnergy);
    scalar Qi = ionizationCrossSections()[maxSpecies].crossSection(eVkinEnergy);
    scalar Qt = Qel+Qex+Qi;

    //Init sigmaTcRMax
    sigmaTcRMax_.primitiveFieldRef() = Qt*cloud.maxwellianMostProbableSpeed
    (
         T,
         cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();
}

template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::elasticCollision(scalar eV, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    CloudType& cloud(this->owner());
    vector preUe, preUn;

    //If the pointer to the neutral species is null we perform a collision with the background gas model
    if(pN != nullptr) {
        preUe = pE->U();
        preUn = pN->U();
        updateVelocity(eV,*pE,*pN);//scatter

        weightCorrection_->correctVelocity(pE,pN,preUe,preUn);//correct weight if necessary

        if(!this->updateNeutralVelocity())//reset if option is selected
            pN->U() = preUn;

    }
    else {//bg collide
        BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

        label celli = pE->cell();
        label typeId = backgroundGas.species();
        vector U = backgroundGas.sampleVelocity(celli);
        vector preU = U;

        updateVelocity(eV,*pE,U,typeId);
        backgroundGas.collisionUpdate(CollisionEvent::Elastic, celli, preU, U);
    }
}

template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::excitationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    CloudType& cloud(this->owner());
    vector preUe, preUn;
    scalar preEV = eV;

    //Excitation, remove the threshold energy
    pE->U() = pE->U() * ::sqrt(1.0 - threshold/eV); //exUe
    eV -= threshold;

    //If the pointer to the neutral species is null we perform a collision with the background gas model
    if(pN != nullptr) {
        preUe = pE->U();
        preUn = pN->U();

        updateVelocity(eV,*pE,*pN);//scatter

        if(isA<NanbuWeightCorrectionModel<CloudType>>(weightCorrection_())) //nanbu weight correction set preUe = pE->U() / ::sqrt(1.0 - threshold/eV)
            preUe /= ::sqrt(1.0 - threshold/preEV);

        weightCorrection_->correctVelocity(pE,pN,preUe,preUn);//correct weight if necessary

        if(!this->updateNeutralVelocity())//reset if option is selected
            pN->U() = preUn;

    }
    else {//bg collide
        BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

        label celli = pE->cell();
        label typeId = backgroundGas.species();
        vector U = backgroundGas.sampleVelocity(celli);
        vector preU = U;

        updateVelocity(eV,*pE,U,typeId);
        backgroundGas.collisionUpdate(CollisionEvent::Exciation, celli, preU, U);
    }
}
template<class CloudType>
bool Foam::ElectronNeutralCollisionModel<CloudType>::ionizationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    bool state = false;//If an Ion was created or not
    vector preUe, preUn;

    scalar massE = cloud.constProps(pE->typeId()).mass();

    //Correct the electrons energy for an ionization event based on the threshold
    scalar E0 = 2.0 - 100.0/(eV+10.0);
    scalar E1 = (eV-threshold)/2.0 - E0;
    scalar alpha = 10.3;
    scalar Ep = E0 + alpha*::tan(rndGen.scalar01()*(::atan(E1/alpha)+::atan(E0/alpha))-::atan(E0/alpha));

    preUe = pE->U();//use "nabbu-ish" weigthing correction, since we can directly use the weight correction model, also needed for new electron
    pE->U() = pE->U()*::sqrt(1.0-(threshold+Ep)/eV);//tilde v as in Nanbu
    eV -= (threshold+Ep);

    //If the pointer to the neutral species is null we perform a collision with the background gas model
    if(pN != nullptr) {

        //Weight fractions used for correction
        scalar wE = pE->nParticle();
        scalar wN = (ionizationNEquivalentParticles_ > 0.0 ? ionizationNEquivalentParticles_ : pN->nParticle());
        scalar wEwN = wE/wN;
        scalar wNwE = 1./wEwN;

        preUn = pN->U();
        updateVelocity(eV,*pE,*pN);//Scatter, updated velocities will be used further down

        scalar rndU = this->owner().rndGen().scalar01();

        //Only create the ion with a propability of wE/wN
        if(rndU <= wEwN) {
            typename CloudType::parcelType* pIon = cloud.addNewParcel(pN->position(),pN->cell(),pN->U(),cloud.ionTypeId(pN->typeId()));//copy already scatted neutral
            pIon->nParticle() = wN;//set weight of new ion to that of the neutral
            state = true;
        }
        if(rndU > wNwE){//prop: rndU < wNwE update electron velocity
            pE->U() = preUe;//so reset velocity
        }

        //Only create the new electron with a propability of wE/wN
        if(rndU <= wEwN) {
            typename CloudType::parcelType* pEnew = cloud.addNewParcel(pE->position(),pE->cell(), \
                            (preUe/mag(preUe))*::sqrt(2.0*Ep*constant::electromagnetic::e.value()/massE),pE->typeId());//tildvp
            pEnew->nParticle() = wN;//set weight of new electron to that of the neutral

            updateVelocity(Ep,*pEnew,preUn,pN->typeId());//Update the velocity if the new electron using energy Ep
        }

        if(!this->updateNeutralVelocity())//in the case we dont delete this is needed
            pN->U() = preUn;

    }
    else {//bg collide
        BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

        label celli = pE->cell();
        label typeId = backgroundGas.species();
        vector U = backgroundGas.sampleVelocity(celli);

        scalar wE = pE->nParticle();
        scalar wN = backgroundGas.nParticle();
        scalar wEwN = wE/wN;
        scalar wNwE = 1./wEwN;

        preUn = U;
        updateVelocity(eV,*pE,U,typeId);
        backgroundGas.collisionUpdate(CollisionEvent::Ionization, celli, preUn, Zero);

        scalar rndU = this->owner().rndGen().scalar01();
        if(rndU <= wEwN) {
            typename CloudType::parcelType* pIon = cloud.addNewParcel(pE->position(),celli,U,cloud.ionTypeId(typeId));//copy already scatted neutral
            pIon->nParticle() = backgroundGas.nParticle();//set weight of new ion to that of the neutral
            state = true;
        }
        if(rndU > wNwE){//prop: rndU < wNwE update electron velocity
            pE->U() = preUe;//so reset velocity
        }
        if(rndU <= wEwN) {
            typename CloudType::parcelType* pEnew = cloud.addNewParcel(pE->position(),celli, \
                            (preUe/mag(preUe))*::sqrt(2.0*Ep*constant::electromagnetic::e.value()/massE),pE->typeId());//tildvp
            pEnew->nParticle() = backgroundGas.nParticle();//set weight of new electron to that of the neutral

            updateVelocity(Ep,*pEnew,preUn,typeId);
        }
    }
    return state;
}

template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::performBackgroundCollisions()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());
    const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label bgTypeId = backgroundGas.species();

    List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());
    label eTypeId = cloud.electronTypeId();

    label elastics = 0;
    label excitations = 0;
    label ionizations = 0;
    label collisionCandidates = 0;

    //update background (may have changed) FIXME: use update flag?
    backgroundCollisionProp_ = 1.0 -exp(-backgroundSigmaTcRMax_*backgroundGas.numberDensity()*cloud.mesh().time().deltaTValue());

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        //Number of electrons
        DynamicList<typename CloudType::parcelType*>& electronParcels(sortedCellOccupancy[celli][eTypeId]);
        label eC = electronParcels.size();

        //Calculate the number of candidates for a collision
        scalar collisions = backgroundCollisionProp_[celli]*eC+backgroundCollisionRemainder_[celli];
        label nCandidates(collisions);
        collisionCandidates+=nCandidates;
        //Save the remainder
        backgroundCollisionRemainder_[celli] = collisions-nCandidates;

        //Go through nCandidates and check if a collision occurs
        for (label i = 0; i < nCandidates; i++)
        {
            //Pick a random electron
            label candidateE = rndGen.sampleAB<label>(0, eC);//FIXME: candidate should only be selected ONCE !!!!!
            typename CloudType::parcelType* electron = electronParcels[candidateE];

            //Calculate the cross section
            scalar massE = cloud.constProps(eTypeId).mass();
            scalar Ue = mag(electron->U());
            scalar gammaE=1.0/::sqrt(1.0-((electron->U()&electron->U())/(constant::universal::c.value()*constant::universal::c.value())));
            scalar eVkinEnergy = (gammaE-1.0)*massE*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());

            scalar sumSigmaT = 0.0;
            scalar sigmaT = backgroundSigmaTcRMax_/Ue;//Ue >> Un => cR == Ue

            scalar Qel = elasticCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qex = excitationCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar exThreshold = excitationCrossSections()[bgTypeId].threshold();
            scalar Qiz = ionizationCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar izThreshold = ionizationCrossSections()[bgTypeId].threshold();

            //Compare against the maximum cross section
            scalar rndU = rndGen.scalar01();
            if(rndU <= (sumSigmaT = Qel)/sigmaT)//elastic collision
            {
                elastics++;
                elasticCollision(eVkinEnergy,electron,nullptr);
                elasticCollisions_[celli] += 1.0;
            }
            else if(eVkinEnergy >= exThreshold && rndU <= (sumSigmaT+=Qex)/sigmaT)//excitation collision
            {
                excitations++;
                excitationCollision(eVkinEnergy,exThreshold,electron,nullptr);
                excitationCollisions_[celli] += 1.0;
            }
            else if(eVkinEnergy >= izThreshold && rndU <= (sumSigmaT+=Qiz)/sigmaT)//ionization collision
            {
                bool ionCreated = ionizationCollision(eVkinEnergy,izThreshold,electron,nullptr);
                if(ionCreated) {
                    ionizations++;
                    ionizationCollisions_[celli] += 1.0;
                }
            }
        }
    }

    //Print some information
    reduce(collisionCandidates, sumOp<label>());
    reduce(elastics, sumOp<label>());
    reduce(excitations, sumOp<label>());
    reduce(ionizations, sumOp<label>());

    scalar dt = mesh.time().deltaTValue();
    scalar beta = dt/(mesh.time().value()-averageStart_);
    ionizationAverage_ = (1.0-beta)*ionizationAverage_ + beta*backgroundGas.nParticle()*ionizations;
    Info << "    Ionization rate (background): " << ionizationAverage_/dt << " particle/s" << endl;

    label collisions =  elastics+excitations+ionizations;
    if (collisionCandidates)
    {
        Info << "    Collisions between electrons and background  = "
             << collisions << nl
             << "    Acceptance rate                              = "
             << scalar(collisions)/scalar(collisionCandidates) << " (candidates = " << collisionCandidates << ")" << nl
             << "    Excitations                                  = "
             << excitations << nl
             << "    Ionizations                                  = "
             << ionizations << nl
             << nl << endl;
    }
    else
        Info<< "    No Collisions between electrons and background" << nl << endl;
}


template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::performCollisions()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());

    List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    label electronTypeId = cloud.electronTypeId();

    scalar deltaT = mesh.time().deltaTValue();

    scalar nParticleMax = 0.0;

    label collisionCandidates = 0;
    label collisions = 0;
    label excitations = 0;
    label ionizations = 0;

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        DynamicList<typename CloudType::parcelType*>& electronParcels(sortedCellOccupancy[celli][electronTypeId]);

        //Get the maximum electron parcel weight in this cell
        forAll(electronParcels,i) {
            typename CloudType::parcelType* p = electronParcels[i];

            scalar nParticle = p->nParticle();
            if(nParticle > nParticleMax)
                nParticleMax = nParticle;
        }
        //Check all neutral species
        forAll(neutralSpecies,i)
        {
            label neutralTypeId = neutralSpecies[i];
            DynamicList<typename CloudType::parcelType*>& cellParcels(sortedCellOccupancy[celli][neutralTypeId]);

            //Get the maximum neutral parcel weight in this cell
            forAll(cellParcels,i) {
                typename CloudType::parcelType* p = cellParcels[i];

                scalar nParticle = p->nParticle();
                if(nParticle > nParticleMax)
                    nParticleMax = nParticle;
            }

            //Do both species exist in this cell?
            label nC = cellParcels.size();
            label eC = electronParcels.size();
            if(nC > 0 && eC > 0)
            {
                scalar sigmaTcRMax = sigmaTcRMax_[celli];

                //Calculate number of candidates
                //FIXME: collisionRemainder_ in the case of multiple neutral species!!!
                scalar selectedPairs =
                        collisionRemainder_[celli]
                        + eC*nC*nParticleMax*sigmaTcRMax*deltaT
                        /mesh.cellVolumes()[celli];

                label nCandidates(selectedPairs);
                collisionRemainder_[celli] = selectedPairs - nCandidates;
                collisionCandidates += nCandidates;

                //Check nCandidates
                for (label c = 0; c < nCandidates; c++)
                {
                    //Pick random parcels
                    label candidateP = rndGen.sampleAB<label>(0, nC);
                    typename CloudType::parcelType* neutral = cellParcels[candidateP];

                    label candidateQ = rndGen.sampleAB<label>(0, eC);
                    typename CloudType::parcelType* electron = electronParcels[candidateQ];

                    //Relative velocity
                    scalar cR = mag(electron->U()-neutral->U());

                    //scalar massN = cloud.constProps(neutralTypeId).mass();
                    scalar massE = cloud.constProps(electronTypeId).mass();
                    //scalar eVkinEnergy = 0.5*massE*cR*cR/(constant::electromagnetic::e.value());

                    //Calc relativistic kin energy
                    scalar gammaE=1.0/::sqrt(1.0-((cR*cR)/(constant::universal::c.value()*constant::universal::c.value())));
                    scalar eVkinEnergy = (gammaE-1.0)*massE*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());

                    //Calculate total cross section
                    scalar Qel = elasticCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qex = excitationCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qi = ionizationCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qt = Qel+Qex+Qi;

                    scalar exThreshold = excitationCrossSections()[neutralTypeId].threshold();
                    scalar izThreshold = ionizationCrossSections()[neutralTypeId].threshold();

                    //Update sigmaTcRMax_ if the current one is large
                    if(Qt*cR > sigmaTcRMax_[celli]) {
                        sigmaTcRMax_[celli] = Qt*cR;
                    }

                    //Check if a collision occurs
                    scalar rndU = rndGen.scalar01();
                    if(rndU <= Qel*cR/sigmaTcRMax)//elastic collision
                    {
                        collisions++;
                        elasticCollision(eVkinEnergy, electron, neutral);
                        elasticCollisions_[celli] += 1.0;
                    }
                    else if(eVkinEnergy >= exThreshold && rndU <= (Qel+Qex)*cR/sigmaTcRMax)//excitation collision
                    {
                        collisions++;
                        excitations++;
                        scalar Eth = excitationCrossSections()[neutralTypeId].threshold();
                        excitationCollision(eVkinEnergy, Eth, electron, neutral);
                        excitationCollisions_[celli] += 1.0;
                    }
                    else if(eVkinEnergy >= izThreshold && rndU <= Qt*cR/sigmaTcRMax)//ionization collision
                    {
                        collisions++;

                        scalar Eth = ionizationCrossSections()[neutralTypeId].threshold();
                        bool ionCreated = ionizationCollision(eVkinEnergy, Eth, electron, neutral);
                        if(ionCreated)
                        {
                            ionizations++;
                            ionizationCollisions_[celli] += 1.0;
                            if(deleteNeutrals())//should neutrals be removed? (can be specified)
                            {
                                //Update the cell parcel list
                                typename DynamicList<typename CloudType::parcelType*>::iterator iter = cellParcels.begin()+candidateP;
                                cellParcels.erase(iter);
                                cloud.deleteParticle(*neutral);
                                nC--;
                                if(nC < 1)
                                    break;
                            }
                        }
                    }
                }
            }

        }
    }

    //Print some information
    reduce(collisions, sumOp<label>());
    reduce(collisionCandidates, sumOp<label>());
    reduce(excitations, sumOp<label>());
    reduce(ionizations, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions between neutral and electron species  = "
                << collisions << nl
                << "    Acceptance rate                     = "
                << scalar(collisions)/scalar(collisionCandidates) << nl
                << "    Excitations                               = "
                << excitations << nl
                << "    Ionizations                               = "
                << ionizations << nl
                << nl << endl;
    }
    else
    {
        Info<< "    No collisions between neutral and electron species" /*(sum of selectedPairs = " << sumSelectedPairs << ")"*/ << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ElectronNeutralCollisionModelNew.C"

// ************************************************************************* //
