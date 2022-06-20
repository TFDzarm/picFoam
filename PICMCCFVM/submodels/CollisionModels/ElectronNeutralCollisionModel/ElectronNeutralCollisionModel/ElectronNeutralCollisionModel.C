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
    backgroundSigmaTcRMax_(
        IOobject
        (
            "picSigmaTcRMax_ElectronBackground",
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
            "backgroundCollisionRemainder_ElectronNeutral",
            owner.mesh().time().timeName(),
            owner.mesh()
        ),
        owner.mesh(),
        dimensionedScalar("backgroundCollisionRemainder_ElectronNeutral", dimless, 0)
   ),
   sigmaTcRMax_(),
   collisionRemainder_(),
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
   ),
   createElectron_(true),
   createIon_(true)
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
    backgroundSigmaTcRMax_(
        IOobject
        (
            "picSigmaTcRMax_ElectronBackground",
            owner.mesh().time().timeName(),
            owner.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        owner.mesh(),
        dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), Zero),
        calculatedFvPatchScalarField::typeName
    ),
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
   sigmaTcRMax_(),
   collisionRemainder_(),
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
           IOobject::READ_IF_PRESENT,
           IOobject::AUTO_WRITE
       ),
       owner.mesh(),
       dimensionedScalar("zero",  dimless, 0.0),
       calculatedFvPatchScalarField::typeName
  ),
  excitationCollisions_
  (
      IOobject
      (
          "excitationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero",  dimless, 0.0),
      calculatedFvPatchScalarField::typeName
  ),
  ionizationCollisions_
  (
      IOobject
      (
          "ionizationElectronCollisions",
          owner.mesh().time().timeName(),
          owner.mesh(),
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ),
      owner.mesh(),
      dimensionedScalar("zero",  dimless, 0.0),
      calculatedFvPatchScalarField::typeName
  ),
  createElectron_(coeffDict_.lookupOrDefault("createElectron",true)),
  createIon_(coeffDict_.lookupOrDefault("createIon",true))
{
    //Create a scalarField per neutral species  for sigmaTcRMax and the collision remainder
    const List<label>& neutralSpecies(owner.neutralSpecies());
    sigmaTcRMax_.setSize(neutralSpecies.size());
    collisionRemainder_.setSize(neutralSpecies.size());

    //Read fields if present or initalize to zero
    forAll(neutralSpecies,i)
    {
        word species = owner.typeIdList()[neutralSpecies[i]];

        Info << "|->    Reading picSigmaTcRMax_Electron" << species << " if present..." << nl;
        sigmaTcRMax_.set(i,
            new volScalarField
            (
                        IOobject
                        (
                            "picSigmaTcRMax_Electron"+species,
                            owner.mesh().time().timeName(),
                            owner.mesh(),
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        owner.mesh(),
                        dimensionedScalar("zero",  dimensionSet(0,3,-1,0,0,0,0), Zero),
                        calculatedFvPatchScalarField::typeName
             )
       );
       Info << "       |= sigmaTcRmax: " << gMax(sigmaTcRMax_[i].primitiveFieldRef()) << endl;

       collisionRemainder_.set(i,
            new volScalarField::Internal
            (
                         IOobject
                         (
                              "collisionRemainder_Electron"+species,
                              owner.mesh().time().timeName(),
                              owner.mesh()
                              ),
                              owner.mesh(),
                              dimensionedScalar("collisionRemainder", dimless, 0.0)
             )
       );
    }
    const BackgroundGasModel<CloudType>& backgroundGas(owner.backgroundGas());

    //Print info about collisions with the the background gas
    if(backgroundGas.active())
    {
        label bgTypeId = backgroundGas.species();

        Info << "|->    BackgroundGasModel " << backgroundGas.type() << " electron collision properties:" << nl
             << "           species: " << owner.typeIdList()[bgTypeId] << nl
             << "           avg number density: " << gAverage(backgroundGas.numberDensity()) << " m^-3" << nl
             << "           avg temperature: " << gAverage(backgroundGas.temperature()) << " K" << nl
             << "           sigmaTcRmax: " << gMax(backgroundSigmaTcRMax_.primitiveFieldRef()) << nl << endl;
    }

    //Initialize the remainder to a random value
    forAll(neutralSpecies,i)
    {
        forAll(collisionRemainder_, celli)
        {
            collisionRemainder_[i][celli] = owner.rndGen().scalar01();
        }
    }
    forAll(backgroundCollisionRemainder_, celli)
    {
        backgroundCollisionRemainder_[celli] = owner.rndGen().scalar01();
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

    Info << "ElectronNeutral CollisionModel:" << endl;

    CloudType& cloud(this->owner());
    const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    const List<label>& neutralSpecies(cloud.neutralSpecies());
    if(neutralSpecies.size() == 0)
    {
        Info << "|->    WARNING ElectronNeutralCollisionModel::initialize() no neutral species defined!" << endl;
        return;
    }

    //Check which mode we use for initialization
    dictionary& crossSectionDict(cloud.picInitialiseDict().subDict("CrossSectionInitialization"));
    word initMode = crossSectionDict.lookup("electronNeutral");
    if(initMode != "FindMaximum" && initMode != "Temperature")
        FatalErrorInFunction << "Invalid initialization mode, valid options are: (FindMaximum, Temperature)" << abort(FatalError);

    bool findMaxium = false;
    if(initMode == "FindMaximum")
        findMaxium = true;

    const typename CloudType::parcelType::constantProperties& cP(cloud.constProps(cloud.electronTypeId()));
    scalar massE = cP.mass();

    //Get the parameter maximumEnergy and energyStep from the initialization dict and search for the maximum cross section
    if(findMaxium)
    {
        //For all particle species...
        dictionary& searchParamDict(crossSectionDict.subDict("electronNeutral_searchParameter"));
        scalar maximumEnergy = readScalar(searchParamDict.lookup("maximumEnergy"));
        scalar energyStep = readScalar(searchParamDict.lookup("energyStep"));

        Info << "|->    Looking up and initializing sigmaTcRMax to the maximum value in range 0.0 to " << maximumEnergy << " eV..." << endl;
        forAll(neutralSpecies,i)
        {
            label typeId = neutralSpecies[i];

            scalar eV = 0.0;
            scalar sigmaTMax = 0.0;
            for(;eV < maximumEnergy; eV+=energyStep)
            {
                scalar Qel = elasticCrossSections()[typeId].crossSection(eV);
                scalar Qex = excitationCrossSections()[typeId].crossSection(eV);
                scalar Qiz = ionizationCrossSections()[typeId].crossSection(eV);

                sigmaTMax = max(sigmaTMax,::sqrt(2.0*constant::electromagnetic::e.value()*eV/massE)*(Qel+Qex+Qiz));
            }
            //Store the found maximum
            sigmaTcRMax_[i].primitiveFieldRef() = sigmaTMax;
        }
        //For the background gas
        const BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());
        if(backgroundGas.active())
        {
            label bgTypeId = backgroundGas.species();

            scalar eV = 0.0;
            scalar sigmaTMax = 0.0;
            for(;eV < maximumEnergy; eV+=energyStep)
            {
                scalar Qel = elasticCrossSections()[bgTypeId].crossSection(eV);
                scalar Qex = excitationCrossSections()[bgTypeId].crossSection(eV);
                scalar Qiz = ionizationCrossSections()[bgTypeId].crossSection(eV);

                sigmaTMax = max(sigmaTMax,::sqrt(2.0*constant::electromagnetic::e.value()*eV/massE)*(Qel+Qex+Qiz));
            }
            //Store the found maximum
            backgroundSigmaTcRMax_.primitiveFieldRef() = sigmaTMax;
        }
    }
    else //Use the specified temperature and initalize the cross section
    {
        //For all particle species...
        Info << "|->    Initializing sigmaTcRMax by specified temperatures..." << endl;
        scalar T = temperatures[cloud.electronTypeId()];
        if(T < VSMALL) {
            Info << "|->    WARNING in initialization of picSigmaTcRMax_Electron: Zero temperature for electron species" << nl << "       |= Using default value of 1 eV..." << endl;
            T = 11604.52;
        }

        forAll(neutralSpecies,i)
        {
            label typeId = neutralSpecies[i];
            scalar eVkinEnergy = T*constant::physicoChemical::k.value()/constant::electromagnetic::e.value();

            //Calculate the total cross section
            scalar Qel = elasticCrossSections()[typeId].crossSection(eVkinEnergy);
            scalar Qex = excitationCrossSections()[typeId].crossSection(eVkinEnergy);
            scalar Qi = ionizationCrossSections()[typeId].crossSection(eVkinEnergy);
            scalar Qt = Qel+Qex+Qi;

            //Init sigmaTcRMax
            sigmaTcRMax_[i].primitiveFieldRef() = Qt*cloud.maxwellianMostProbableSpeed
            (
                 T,
                 massE
            );

            sigmaTcRMax_[i].correctBoundaryConditions();

        }

        //For the background gas
        if(backgroundGas.active())
        {
            label bgSpecies = backgroundGas.species();
            scalar eVkinEnergy = T*constant::physicoChemical::k.value()/constant::electromagnetic::e.value();

            //Calculate the total cross section
            scalar Qel = elasticCrossSections()[bgSpecies].crossSection(eVkinEnergy);
            scalar Qex = excitationCrossSections()[bgSpecies].crossSection(eVkinEnergy);
            scalar Qi = ionizationCrossSections()[bgSpecies].crossSection(eVkinEnergy);
            scalar Qt = Qel+Qex+Qi;

            //Init sigmaTcRMax
            backgroundSigmaTcRMax_.primitiveFieldRef() = Qt*cloud.maxwellianMostProbableSpeed
            (
                 T,
                 massE
            );
            backgroundSigmaTcRMax_.correctBoundaryConditions();
        }
    }
}

//Elastic Collision of particles
template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::elasticCollision(scalar eV, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    vector preUe = pE->U();
    vector preUn = pN->U();
    updateVelocity(eV,*pE,*pN);//scatter

    weightCorrection_->correctVelocity(pE,pN,preUe,preUn);//correct weight if necessary

    if(!this->updateNeutralVelocity())//reset if option is selected
        pN->U() = preUn;

}


//Elastic Collision with the background gas
template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::elasticCollision(scalar eV, typename CloudType::parcelType* pE, vector& Un)
{
    CloudType& cloud(this->owner());

    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label celli = pE->cell();
    label typeId = backgroundGas.species();
    vector preU = Un;

    updateVelocity(eV,*pE,Un,typeId);//scatter
    backgroundGas.collisionUpdate(CollisionEvent::Elastic, celli, preU, Un);//notify the model of this collision event
}


//Excitation Collision of particles
template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::excitationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    vector preUn = pN->U();
    scalar preEV = eV;

    //Excitation, remove the threshold energy
    pE->U() = pE->U() * ::sqrt(1.0 - threshold/eV);
    eV -= threshold;

    vector preUe = pE->U();

    updateVelocity(eV,*pE,*pN);//scatter

    if(isA<NanbuWeightCorrectionModel<CloudType>>(weightCorrection_())) //nanbu weight correction set preUe = pE->U() / ::sqrt(1.0 - threshold/eV)
        preUe /= ::sqrt(1.0 - threshold/preEV);

    weightCorrection_->correctVelocity(pE,pN,preUe,preUn);//correct weight if necessary

    if(!this->updateNeutralVelocity())//reset if option is selected
        pN->U() = preUn;

}

//Excitation Collision with the background gas
template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::excitationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, vector& Un)
{
    CloudType& cloud(this->owner());

    //Excitation, remove the threshold energy
    pE->U() = pE->U() * ::sqrt(1.0 - threshold/eV);
    eV -= threshold;


    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label celli = pE->cell();
    label typeId = backgroundGas.species();
    vector preU = Un;

    updateVelocity(eV,*pE,Un,typeId);//scatter
    backgroundGas.collisionUpdate(CollisionEvent::Exciation, celli, preU, Un);//notify the model of this collision event
}

//Ionization Collision of particles
template<class CloudType>
bool Foam::ElectronNeutralCollisionModel<CloudType>::ionizationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, typename CloudType::parcelType* pN)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    bool state = false;//If an Ion was created or not
    vector preUn = pN->U();

    scalar massE = cloud.constProps(pE->typeId()).mass();

    //Correct the electrons energy for an ionization event based on the threshold
    scalar E0 = 2.0 - 100.0/(eV+10.0);
    scalar E1 = (eV-threshold)/2.0 - E0;
    scalar alpha = 10.3;
    scalar Ep = E0 + alpha*::tan(rndGen.scalar01()*(::atan(E1/alpha)+::atan(E0/alpha))-::atan(E0/alpha));

    vector preUe = pE->U();//use "nabbu-ish" weigthing correction, also needed for new electron
    pE->U() = pE->U()*::sqrt(1.0-(threshold+Ep)/eV);//tildv as in Nanbu
    eV -= (threshold+Ep);

    //Weight fractions used for correction
    scalar wE = pE->nParticle();
    scalar wN = (ionizationNEquivalentParticles_ > 0.0 ? ionizationNEquivalentParticles_ : pN->nParticle());
    scalar wEwN = wE/wN;
    scalar wNwE = 1./wEwN;

    updateVelocity(eV,*pE,*pN);//Scatter, updated velocities will be used further down

        //Alternative: Always create the ion/electron pair and reduce the neutral mass (may fill up the domain with lots of particles depending on weight difference)
        /*
        typename CloudType::parcelType* pIon = cloud.addNewParcel(pN->position(),pN->cell(),pN->U(),cloud.ionTypeId(pN->typeId()));
        pIon->nParticle() = wE;//Weight of the new ion is weight of the electron
                               //What should happen if wN < wE ... set wI und wEnew to wN... but than the post collision velocities need adjustment...
        if(this->deleteNeutrals())//Only reduce the mass if we also want the delete the neutral eventually
        {
            pN->nParticle() -= wE;//Reduce weight of neutral
            //At this point an adjustment for the lost momentum in the neutral is missing...
        }
        state = true;

        typename CloudType::parcelType* pEnew = cloud.addNewParcel(pE->position(),pE->cell(), \
                        (preUe/mag(preUe))*::sqrt(2.0*Ep*constant::electromagnetic::e.value()/massE),pE->typeId());//tildvp
        pEnew->nParticle() = wE;//set weight of new electron to that of the electron
        updateVelocity(Ep,*pEnew,preUn,pN->typeId());
        */
        //END Alternative


    //Only create the ion with a propability of wE/wN and if the option createIonElectronPair() is set to true
    scalar rndU = cloud.rndGen().scalar01();
    if(rndU <= wEwN && createIon()) {
        typename CloudType::parcelType* pIon = cloud.addNewParcel(pN->position(),pN->cell(),pN->U(),cloud.ionTypeId(pN->typeId()));//copy already scatted neutral
        pIon->nParticle() = wN;//set weight of new ion to that of the neutral
        state = true;
    }

    if(rndU > wNwE){//prop: rndU < wNwE update electron velocity
        pE->U() = preUe;//so reset velocity
    }

    //Only create the new electron with a propability of wE/wN and if the option createIonElectronPair() is set to true
    if(rndU <= wEwN && createElectron()) {
        typename CloudType::parcelType* pEnew = cloud.addNewParcel(pE->position(),pE->cell(), \
                        (preUe/mag(preUe))*::sqrt(2.0*Ep*constant::electromagnetic::e.value()/massE),pE->typeId());//tildvp
        pEnew->nParticle() = wN;//set weight of new electron to that of the neutral
        updateVelocity(Ep,*pEnew,preUn,pN->typeId());//Update the velocity if the new electron using energy Ep
    }

    if(!this->updateNeutralVelocity())//in the case we dont delete this is needed
        pN->U() = preUn;

    return state;
}


//Ionization Collision with the background gas
template<class CloudType>
bool Foam::ElectronNeutralCollisionModel<CloudType>::ionizationCollision(scalar eV, scalar threshold, typename CloudType::parcelType* pE, vector& Un)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());

    bool state = false;//ion created or not
    vector preUn = Un;

    scalar massE = cloud.constProps(pE->typeId()).mass();

    //Correct the electrons energy for an ionization event based on the threshold
    scalar E0 = 2.0 - 100.0/(eV+10.0);
    scalar E1 = (eV-threshold)/2.0 - E0;
    scalar alpha = 10.3;
    scalar Ep = E0 + alpha*::tan(rndGen.scalar01()*(::atan(E1/alpha)+::atan(E0/alpha))-::atan(E0/alpha));

    vector preUe = pE->U();//use "nabbu-ish" weigthing correction, also needed for new electron
    pE->U() = pE->U()*::sqrt(1.0-(threshold+Ep)/eV);//tildv
    eV -= (threshold+Ep);


    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label celli = pE->cell();
    label typeId = backgroundGas.species();

    //Weight fractions used for correction
    scalar wE = pE->nParticle();
    scalar wN = backgroundGas.nParticle();
    scalar wEwN = wE/wN;
    scalar wNwE = 1./wEwN;

    updateVelocity(eV,*pE,Un,typeId);//Scatter, updated velocities will be used further down
    backgroundGas.collisionUpdate(CollisionEvent::Ionization, celli, preUn, Un);//notify the model of the collision event

    //Only create the ion with a propability of wE/wN and if the option createIonElectronPair() is set to true
    scalar rndU = cloud.rndGen().scalar01();
    if(rndU <= wEwN && createIon()) {
        typename CloudType::parcelType* pIon = cloud.addNewParcel(pE->position(),celli,Un,cloud.ionTypeId(typeId));//copy already scatted neutral
        pIon->nParticle() = backgroundGas.nParticle();//set weight of new ion to that of the neutral
        state = true;
    }

    if(rndU > wNwE){//prop: rndU < wNwE update electron velocity
        pE->U() = preUe;//so reset velocity
    }

    //Only create the new electron with a propability of wE/wN and if the option createIonElectronPair() is set to true
    if(rndU <= wEwN && createElectron()) {
        typename CloudType::parcelType* pEnew = cloud.addNewParcel(pE->position(),celli, \
                        (preUe/mag(preUe))*::sqrt(2.0*Ep*constant::electromagnetic::e.value()/massE),pE->typeId());//tildvp
        pEnew->nParticle() = backgroundGas.nParticle();//set weight of new electron to that of the neutral
        updateVelocity(Ep,*pEnew,preUn,typeId);
    }
    return state;
}

template<class CloudType>
void Foam::ElectronNeutralCollisionModel<CloudType>::performBackgroundCollisions()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());
    BackgroundGasModel<CloudType>& backgroundGas(cloud.backgroundGas());

    label bgTypeId = backgroundGas.species();

    List<List<DynamicList<typename CloudType::parcelType*>>>& sortedCellOccupancy(cloud.sortedCellOccupancy());
    label eTypeId = cloud.electronTypeId();

    label elastics = 0;
    label excitations = 0;
    label ionizations = 0;
    label collisionCandidates = 0;

    scalar dt = mesh.time().deltaTValue();
    scalar massE = cloud.constProps(eTypeId).mass();
    const scalarField& bgDensity(backgroundGas.numberDensity());

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        //Number of electrons
        DynamicList<typename CloudType::parcelType*>& electronParcels(sortedCellOccupancy[celli][eTypeId]);
        label eC = electronParcels.size();

        scalar sigmaTcRMax = backgroundSigmaTcRMax_[celli];

        //Calculate the number of candidates for a collision
        scalar collisions = backgroundCollisionRemainder_[celli] + eC*bgDensity[celli]*sigmaTcRMax*dt;

        //Proposed in Nanbu 2000, used in xpdp1, not in oopd1...
        //scalar collisions = backgroundCollisionRemainder_[celli] + (1.0 - exp(-sigmaTcRMax*bgDensity[celli]*dt))*eC;

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

            //Sample velocity / Alternative: most probable speed = sqrt(2.0*constant::physicoChemical::k.value()*bgTemperature[celli]/massN)
            vector Un = backgroundGas.sampleVelocity(celli);
            scalar cR = mag(electron->U()-Un);
            scalar gamma=1.0/::sqrt(1.0-((cR*cR)/(constant::universal::c.value()*constant::universal::c.value())));

            //Calc relativistic kin energy
            scalar eVkinEnergy = (gamma-1.0)*massE*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());
            //scalar eVkinEnergy = 0.5*massE*(cR*cR)/constant::electromagnetic::e.value();

            scalar sigmaT = sigmaTcRMax/cR;

            //Calculate the cross section
            scalar Qel = elasticCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qex = excitationCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qiz = ionizationCrossSections()[bgTypeId].crossSection(eVkinEnergy);
            scalar Qt = Qel+Qex+Qiz;

            scalar exThreshold = excitationCrossSections()[bgTypeId].threshold();
            scalar izThreshold = ionizationCrossSections()[bgTypeId].threshold();


            //Update backgroundSigmaTcRMax_... this needs to be done because of fixed/constant cross sections
            if(Qt*cR > backgroundSigmaTcRMax_[celli]) {
                backgroundSigmaTcRMax_[celli] = Qt*cR;
            }

            //Compare against the maximum cross section
            scalar rndU = rndGen.scalar01();
            if(rndU <= Qel/sigmaT) //Qel*cR/sigmaTcRMax
            {
                elastics++;
                elasticCollision(eVkinEnergy,electron,Un);
                elasticCollisions_[celli] += 1.0;
            }
            else if(eVkinEnergy >= exThreshold && rndU <= (Qel+Qex)/sigmaT)//(Qel+Qex)*cR/sigmaTcRMax
            {
                excitations++;
                excitationCollision(eVkinEnergy,exThreshold,electron,Un);
                excitationCollisions_[celli] += 1.0;
            }
            else if(eVkinEnergy >= izThreshold && rndU <= Qt/sigmaT)//Qt*cR/sigmaTcRMax
            {
                bool ionCreated = ionizationCollision(eVkinEnergy,izThreshold,electron,Un);
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

    label collisionCandidates = 0;
    label collisions = 0;
    label excitations = 0;
    label ionizations = 0;

    // Go through all cells
    forAll(sortedCellOccupancy, celli)
    {
        //Get the maximum electron parcel weight in this cell
        scalar nParticleMax = 0.0;
        DynamicList<typename CloudType::parcelType*>& electronParcels(sortedCellOccupancy[celli][electronTypeId]);

        forAll(electronParcels,i) {
            typename CloudType::parcelType* p = electronParcels[i];

            scalar nParticle = p->nParticle();
            if(nParticle > nParticleMax)
                nParticleMax = nParticle;
        }

        //Check all neutral species
        forAll(neutralSpecies,i)
        {
            scalar sigmaTcRMax = sigmaTcRMax_[i][celli];
            if(sigmaTcRMax <= 0.0)
                continue;

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
                //Calculate number of candidates
                scalar selectedPairs =
                        collisionRemainder_[i][celli]
                        + eC*nC*nParticleMax*sigmaTcRMax*deltaT/mesh.cellVolumes()[celli];

                label nCandidates(selectedPairs);
                collisionRemainder_[i][celli] = selectedPairs - nCandidates;
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
                    scalar massE = cloud.constProps(electronTypeId).mass();

                    //Calc relativistic kin energy
                    scalar gammaE=1.0/::sqrt(1.0-((cR*cR)/(constant::universal::c.value()*constant::universal::c.value())));
                    scalar eVkinEnergy = (gammaE-1.0)*massE*constant::universal::c.value()*constant::universal::c.value()/(constant::electromagnetic::e.value());
                    //scalar eVkinEnergy = 0.5*massE*cR*cR/(constant::electromagnetic::e.value());

                    //Calculate total cross section
                    scalar Qel = elasticCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qex = excitationCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qi = ionizationCrossSections()[neutralTypeId].crossSection(eVkinEnergy);
                    scalar Qt = Qel+Qex+Qi;

                    scalar exThreshold = excitationCrossSections()[neutralTypeId].threshold();
                    scalar izThreshold = ionizationCrossSections()[neutralTypeId].threshold();

                    //Update sigmaTcRMax_ if the current one is large
                    if(Qt*cR > sigmaTcRMax_[i][celli]) {
                        sigmaTcRMax_[i][celli] = Qt*cR;
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

                        bool ionCreated = ionizationCollision(eVkinEnergy, izThreshold, electron, neutral);
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
                            /*if(deleteNeutrals() && neutral->nParticle() <= 0.0) //Used for the alternative approch of allways creating ions/electrons, here we need to delete the neutral if the weight is zero
                            {
                                typename DynamicList<typename CloudType::parcelType*>::iterator iter = cellParcels.begin()+candidateP;
                                cellParcels.erase(iter);
                                cloud.deleteParticle(*neutral);
                                nC--;
                                if(nC < 1)
                                    break;
                            }*/
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

    forAll(neutralSpecies,i)
        sigmaTcRMax_[i].correctBoundaryConditions();

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
