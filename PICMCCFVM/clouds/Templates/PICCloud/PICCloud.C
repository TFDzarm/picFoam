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

#include "PICCloud.H"
#include "BinaryCollisionModel.H"
#include "WallReflectionModel.H"
#include "BoundaryModel.H"
#include "MaxwellSolver.H"
#include "ParticlePusher.H"
#include "CoulombCollisionModel.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "polyMeshTetDecomposition.H"
#include "BoundaryEvent.H"
#include "BoundaryEventModelList.H"
#include "ElectronNeutralCollisionModel.H"
#include "IonNeutralCollisionModel.H"
#include "BoundaryModelList.H"
#include "DiagnosticsList.H"
#include "BackgroundGasModel.H"
#include "ParticleMerging.H"
#include "ChargeDistribution.H"
#include "meshTools.H"

using namespace Foam::constant;
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
Foam::PICCloud<ParcelType>::buildConstProps() is called in the PICCloud constructor.

This reads the "moleculeProperties" entry in the constant/picProperties file and initializes the list containing constant properties for each species.
Additionally, lists linking species together e.g. a neutal species with its ion (used in the collision procedures) are initialized.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    //Read the sub-dictionary entry "moleculeProperties" from constant/picProperties
    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    //Print species info, check for equal particle weight
    scalar eqWeight = 0.0;
    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] =
        typename ParcelType::constantProperties(molDict);

        Info<< "    " << id << nl
        << "    |_ mass:          " << constProps_[i].mass() << " [kg]" << nl
        << "    |_ d:             " << constProps_[i].d() << " [m]" << nl
        << "    |_ charge:        " << constProps_[i].charge() << " [C]" << nl
        << "    |_ weight:        " << constProps_[i].nParticle() << " [-]" << nl
        << "    |_ solveMovement: " << (constProps_[i].solveMovement() ? "yes" : "no")
        << (i == typeIdList_.size()-1 ? "\n" : "") << endl;

        if(i == 0)//fist value
            eqWeight = constProps_[i].nParticle();
        else if(constProps_[i].nParticle() != eqWeight)
            eqWeight = 0.0;
    }
    //If all species have the same weight this varibales is set this weight, else it is zero. This used to check if some of collision algorithm that are only applicable to equal weighted particle species can be used.
    nParticleEqWeight_ = eqWeight;

    //Read species relations used in CollisionModels e.g. for ionization
    dictionary relations
    (
        moleculeProperties.subDict("SpeciesRelations")
    );

    word electronType = relations.lookup("electronTypeId");
    electronTypeId_ = findIndex(typeIdList_,electronType);

    //Relations of species: ions -> neutral, neutral -> ions, not existing == -1
    neutralTypeList_.setSize(typeIdList_.size(),-1);
    ionTypeList_.setSize(typeIdList_.size(),-1);

    forAllConstIter(IDLList<entry>, relations, iter)
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();
            word nSType = subDict.lookup("neutralTypeId");
            label nS = findIndex(typeIdList_,nSType);

            word iSType = subDict.lookup("ionTypeId");
            label iS = findIndex(typeIdList_,iSType);

            if(nS != -1)
            {
                neutralTypeList_[nS] = nS;
                ionTypeList_[nS] = iS;
            }

            if(iS != -1)
            {
                neutralTypeList_[iS] = nS;
                ionTypeList_[iS] = iS;
            }
        }
    }

    //Lists for different species groups...(ions, neutrals, charged)
    DynamicList<label> cS;
    DynamicList<label> nS;
    DynamicList<label> iS;
    forAll(typeIdList_,typeId)
    {
        if(constProps(typeId).charge() != 0.0)
        {
            cS.append(typeId);
            if(typeId != electronTypeId_)
                iS.append(typeId);
        }
        else if(constProps(typeId).charge() == 0.0)
            nS.append(typeId);
    }

    chargedSpecies_.transfer(cS);
    neutralSpecies_.transfer(nS);
    ionSpecies_.transfer(iS);

    //Check if we want to calculate fields like number density separately for these given species
    wordList fieldCalc = particleProperties_.subDict("SolverSettings").lookup("fieldCalculation");
    label index = 0;
    forAll(fieldCalc,i)
    {
        word species = fieldCalc[i];
        label typeId = findIndex(typeIdList_,species);
        if(typeId < 0)
            FatalErrorInFunction << "fieldCalculation: no species " << species << " defined" << endl;

        fieldCalculation_[typeId] = index++;
    }

    if(index > 0)
    {
        //Initialize the fields for each species definied in "fieldCalculation"
        N_.setSize(index);
        momentumSpecies_.setSize(index);
        linearKESpecies_.setSize(index);
        rhoNSpecies_.setSize(index);
        rhoMSpecies_.setSize(index);
        rhoChargeSpecies_.setSize(index);
        jSpecies_.setSize(index);

        for(label fieldId = 0; fieldId < index; fieldId++)
        {
            label typeId = findIndex(fieldCalculation_,fieldId);
            word species = typeIdList_[typeId];

            N_.set(fieldId,
                new volScalarField
                (
                            IOobject
                            (
                                "N:"+species,
                                mesh_.time().timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh_,
                            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0, 0,0), Zero),
                            calculatedFvPatchScalarField::typeName
                 )
           );

           momentumSpecies_.set(fieldId,
                         new volVectorField(
                            IOobject
                            (
                                "momentum:"+species,
                                mesh_.time().timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh_,
                            dimensionedVector
                            (
                                "zero",
                                dimensionSet(1, -2, -1, 0, 0),
                                Zero
                            ),
                            calculatedFvPatchVectorField::typeName
                         )
                       );

           linearKESpecies_.set(fieldId,
             new volScalarField
             (
               IOobject
               (
                   "linearKE:"+species,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh_,
               dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0),
               calculatedFvPatchVectorField::typeName
           ));
           rhoNSpecies_.set(fieldId,
             new volScalarField
             (
               IOobject
               (
                   "rhoN:"+species,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh_,
               dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0),
               calculatedFvPatchVectorField::typeName
           ));
           rhoMSpecies_.set(fieldId,
             new volScalarField
             (
               IOobject
               (
                   "rhoM:"+species,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh_,
               dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), 0.0),
               calculatedFvPatchVectorField::typeName
           ));

           rhoChargeSpecies_.set(fieldId,
           new volScalarField
           (
               IOobject
               (
                   "rhoCharge:"+species,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh_,
               dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0,1,0), 0.0),
               zeroGradientFvPatchScalarField::typeName
           ));

           jSpecies_.set(fieldId,
           new volVectorField
           (
               IOobject
               (
                   "j:"+species,
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh_,
               dimensionedVector("zero",  dimensionSet(0, -2, 0, 0, 0,1,0), Zero),
               calculatedFvPatchScalarField::typeName
           ));

        }
    }
}


/*
Foam::PICCloud<ParcelType>::buildCellOccupancy() is called in the evolve function and the constructor.

This updates the field sortedCellOccupancy_ which containes the pointers to each parcel contained within specific cells.
Additionally the density fields of species defined by the fieldCalculation are updated.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::buildCellOccupancy()
{
    //Clear species fieldCalculation density fields
    forAll(N_,i)
    {
        N_[i] = dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0), 0.0);
    }

    //Clear species sorted cell occupancy
    forAll(sortedCellOccupancy_, cO)
    {
        forAll(typeIdList_, typeId)
        {
            sortedCellOccupancy_[cO][typeId].clear();
        }
    }

    //Calculate cell occupancy and fieldCalculation density
    forAllIter(typename PICCloud<ParcelType>, *this, iter)
    {
        sortedCellOccupancy_[iter().cell()][iter().typeId()].append(&iter());

        label fId = fieldCalculation_[iter().typeId()];
        if(fId >= 0) {
            N_[fId][iter().cell()] += iter().nParticle();
        }
    }
}


/*
Foam::PICCloud<ParcelType>::setupParticleProperties() is called in the constructor.

This calculates the charge of each particle once at he beginning of the simulation since we only save the ionization state not the actual charge value of each parcel.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::setupParticleProperties()
{
    //Setup particle charge: q=q_init*Z
    forAllIter(typename PICCloud<ParcelType>, *this, iter)
    {
        iter().charge() = constProps(iter().typeId()).charge()*iter().chargeModifier();
    }
}

/*
Foam::PICCloud<ParcelType>::setModels() is called in the constructor.

This initializes all submodels and verifies if the simulation uses equal weighted species or not.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::setModels()
{
    Info << "Initializing submodels" << nl << "----------------------" << endl;
    //Pick the charge distribution model: Charge interpolation from parcel to the cells. 
    chargeDensity_.reset(
        ChargeDistribution<PICCloud<ParcelType>>::New(
            "rhoCharge",
            particleProperties_.subDict("SolverSettings"),
            *this
        ).ptr()
    );
    //Pick the interpolation model from the electric field to the parcels.
    eField_.reset(
            FieldWeigthing::New(
            "E",
            particleProperties_.subDict("SolverSettings"),
            mesh_
        ).ptr()
    );
    //Pick the background collision model.
    //Load this one before collision models!
    backgroundGas_.reset(
        BackgroundGasModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("CollisionModels"),
            *this
        ).ptr()
    );
    //Pick the merging algortihm.
    particleMerging_.reset(
        ParticleMerging<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        ).ptr()
    );

    //Check if the simulation uses equal weight.
    if(particleMerging_->active())
        nParticleEqWeight_ = 0.0;//if we merge particles there cannot be equal weigths...
    if(!particleMerging_->active() && backgroundGas_->active() && backgroundGas_->nParticle() != nParticleEqWeight_)//Background model with different weight
        nParticleEqWeight_ = 0.0;

    //Check before collision models are initilized.... !!!
    //FIXME: Does not check for emitter in BoundaryModels!!!
    bool check = nParticleEqWeight_ != 0.0;
    if(check) {
        Info << "+ Checking particle weight..." << nl;
        forAllIter(typename PICCloud<ParcelType>, *this, iter)
        {
            ParcelType& p = iter();
            if(p.nParticle() != nParticleEqWeight_)
            {
                nParticleEqWeight_ = 0.0;
                check = false;
                break;
            }
        }
        if(check)
            Info << "|->    Simulation uses equal weighted particles" << nl << endl;
    }

    if(!check)
        Info << "|->    Particles have different weighting factors" << nl << endl;


    //Pick the binary collision model (neutal-neutal and neutal-ion).
    binaryCollisionModel_.reset(
        BinaryCollisionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("CollisionModels"),
            *this
        ).ptr()
    );
    //Pick the coulomb collision model (charged-charged).
    coulombCollisionModel_.reset(
        CoulombCollisionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("CollisionModels"),
            *this
        ).ptr()
    );
    //Pick the electron-neutral collision model.
    electronNeutralCollisionModel_.reset(
        ElectronNeutralCollisionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("CollisionModels"),
            *this
        ).ptr()
    );
    //Pick the ion-neutral collision model.
    ionNeutralCollisionModel_.reset(
        IonNeutralCollisionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("CollisionModels"),
            *this
            ).ptr()
    );
    //Pick the solver for the electric field.
    maxwellSolver_.reset(
        MaxwellSolver<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        ).ptr()
    );
    //Pick the particle pusher algortihm.
    particlePusher_.reset(
        ParticlePusher<PICCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        ).ptr()
    );

    //Pick the wall interaction model.
    wallReflectionModel_.reset(
        WallReflectionModel<PICCloud<ParcelType>>::New
        (
            particleProperties_.subDict("BoundaryModels"),
            *this
        ).ptr()
    );
    
    //Pick boundary models for patches e.g. circuit models.
    boundaryModels_.setupModels(
        particleProperties_.subDict("BoundaryModels").subDict("PatchBoundaryModels"),
        *this
    );
    //Boundary events, mainly for diagnotic purpose.
    boundaryEvents_.setupModels(
        particleProperties_.subDict("BoundaryModels").subOrEmptyDict("PatchEventModels"),
        *this
    );
    //Pick diagnostics models.
    parcelDiagnostics_.setupModels(
        particleProperties_.subDict("Diagnostics"),
        *this
    );

}

/*
Foam::PICCloud<ParcelType>::collisions() is called in the evolve function.

Call all collision submodels.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::collisions()
{
    if(coulombCollision().active())
        coulombCollision().handleCollisions();

    if (electronNeutralCollision().active())
        electronNeutralCollision().handleCollisions();

    if (ionNeutralCollision().active())
        ionNeutralCollision().handleCollisions();

    if (binaryCollision().active())
        binaryCollision().handleCollisions();
}

/*
Foam::PICCloud<ParcelType>::resetFields() is called in the evolve function.

Reset fields calculated every time step from the parcel distribution.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::resetFields()
{
    q_ = dimensionedScalar( dimensionSet(1, 0, -3, 0, 0), 0);

    fD_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -1, -2, 0, 0),
        Zero
    );

    rhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), vSmall);
    rhoM_ =  dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), vSmall);
    picRhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0);
    linearKE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    momentum_ = dimensionedVector
    (
         "zero",
         dimensionSet(1, -2, -1, 0, 0),
         Zero
    );

    forAll(momentumSpecies_,i)
    {
        rhoNSpecies_[i] = dimensionedScalar("zero",dimensionSet(0, -3, 0, 0, 0),0.0);
        rhoMSpecies_[i] = dimensionedScalar("zero",dimensionSet(1, -3, 0, 0, 0),0.0);
        momentumSpecies_[i] = dimensionedVector("zero",dimensionSet(1, -2, -1, 0, 0),Zero);
        linearKESpecies_[i] = dimensionedScalar("zero",dimensionSet(1, -1, -2, 0, 0),0.0);
    }
    //Do not reset j_ and rhoCharge_!!! This is done in calculateFields()
}

/*
Foam::PICCloud<ParcelType>::calculateFields() is called in the evolve function.

Calculated fields every time step from the parcel distribution. Some field are resetted here, which are needed in between calls to resetFields() and this functions.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::calculateFields()
{
    //Reset j_ and the charge density
    ChargeDistribution<PICCloud<ParcelType>>& cD = chargeDensity_();
    cD.reset();

    j_ = dimensionedVector
    (
        "zero",
        dimensionSet(0, -2, 0, 0, 0, 1, 0),
        Zero
    );

    //Reset fieldCalculation fields
    forAll(rhoChargeSpecies_,i)
    {
        rhoChargeSpecies_[i] = dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0,1,0), 0.0);
        jSpecies_[i] = dimensionedVector("zero",dimensionSet(0, -2, 0, 0, 0, 1, 0),Zero);
    }


    vectorField& j = j_.primitiveFieldRef();

    scalarField& rhoN = rhoN_.primitiveFieldRef();
    scalarField& rhoM = rhoM_.primitiveFieldRef();
    scalarField& picRhoN = picRhoN_.primitiveFieldRef();
    scalarField& linearKE = linearKE_.primitiveFieldRef();
    vectorField& momentum = momentum_.primitiveFieldRef();

    //Go through the linked list and calculate the properties...
    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label celli = p.cell();
        const label fId = fieldCalculation_[p.typeId()];

        //Add particle to charge density model
        cD.add(p);

        //Number of molecules
        rhoN[celli] += p.nParticle();
        //Mass of the parcels
        rhoM[celli] += p.nParticle()*constProps(p.typeId()).mass();
        //Number of parcel
        picRhoN[celli]++;
        //Linear kinetic energy
        linearKE[celli] += 0.5*p.nParticle()*constProps(p.typeId()).mass()*(p.U() & p.U());
        //Momentum
        momentum[celli] += p.nParticle()*constProps(p.typeId()).mass()*p.U();

        //Calculated species specific fields.
        if(fId >= 0) {
            jSpecies_[fId][celli] += p.nParticle()*p.charge()*p.U();
            rhoChargeSpecies_[fId][celli] += p.nParticle()*p.charge();

            rhoNSpecies_[fId][celli] += p.nParticle();
            rhoMSpecies_[fId][celli] += p.nParticle()*constProps(p.typeId()).mass();
            momentumSpecies_[fId][celli] += p.nParticle()*constProps(p.typeId()).mass()*p.U();
            linearKESpecies_[fId][celli] += 0.5*p.nParticle()*constProps(p.typeId()).mass()*(p.U() & p.U());
        }

        j[celli] += p.nParticle()*p.charge()*p.U();
    }

    //Divide by the cell volume field.
    rhoN /= mesh().cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM /= mesh().cellVolumes();
    rhoM_.correctBoundaryConditions();

    picRhoN_.correctBoundaryConditions();

    linearKE /= mesh().cellVolumes();
    linearKE_.correctBoundaryConditions();

    momentum /= mesh().cellVolumes();
    momentum_.correctBoundaryConditions();

    forAll(momentumSpecies_,i)
    {
        rhoChargeSpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        rhoChargeSpecies_[i].correctBoundaryConditions();

        jSpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        jSpecies_[i].correctBoundaryConditions();

        rhoNSpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        rhoNSpecies_[i].correctBoundaryConditions();

        rhoMSpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        rhoMSpecies_[i].correctBoundaryConditions();

        momentumSpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        momentumSpecies_[i].correctBoundaryConditions();

        linearKESpecies_[i].primitiveFieldRef() /= mesh().cellVolumes();
        linearKESpecies_[i].correctBoundaryConditions();
    }

    //Finalize charge density model update and add space charge density
    cD.update();
    cD.field() += spaceChargeDensity_;
    cD.field().correctBoundaryConditions();

    j /= mesh().cellVolumes();
    j_.correctBoundaryConditions();

}



// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

/*
Foam::PICCloud<ParcelType>::addNewParcel

This function is used to add a parcel to the cloud.
*/
template<class ParcelType>
ParcelType* Foam::PICCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const label celli,
    const vector& U,
    const label typeId
)
{
    ParcelType* p = new ParcelType(mesh_, position, celli, U, typeId);

    //Setup charge and particle weight
    p->charge() = constProps(typeId).charge();
    p->nParticle() = constProps(typeId).nParticle();

    this->addParticle(p);//PICCloud ist a double linked list, add the particle
    return p;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
Foam::PICCloud<ParcelType>::PICCloud

The constructor: Setup all members of the class.
*/
template<class ParcelType>
Foam::PICCloud<ParcelType>::PICCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.subDict("SolverSettings").lookup("typeIdList")),
    nParticleEqWeight_(0.0),
    chargedSpecies_(),
    neutralSpecies_(),
    ionSpecies_(),
    neutralTypeList_(),
    ionTypeList_(),
    electronTypeId_(-1),
    fieldCalculation_(typeIdList_.size(),-1),
    sortedCellOccupancy_(mesh_.nCells()),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    fD_
    (
        IOobject
        (
            "fD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoNSpecies_(),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoMSpecies_(),
    rhoM_
    (
        IOobject
        (
            "rhoM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    picRhoN_
    (
        IOobject
        (
            "picRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    linearKESpecies_(),
    linearKE_
    (
        IOobject
        (
            "linearKE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    momentumSpecies_(),
    momentum_
    (
        IOobject
        (
            "momentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoChargeSpecies_(),
    chargeDensity_(),
    N_(),
    eField_(),
    phiE_
        (
            IOobject
            (
                "phiE",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        ),
    jSpecies_(),
    j_
    (
        IOobject
        (
            "j",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero",  dimensionSet(0, -2, 0, 0, 0,1,0), Zero),
        calculatedFvPatchScalarField::typeName
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    jE_
    (
        IOobject
        (
            "jE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -3, 0, 0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    spaceChargeDensity_
    (
        IOobject
        (
            "spaceChargeDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0,1,0), Zero)
    ),
    cellLengthScale_(mesh_.nCells(),0.0),
    avgcellLengthScale_(0.0),
    printVelocityWarning_(true),
    syncVelocityAtBoundary_(particleProperties_.subDict("SolverSettings").lookupOrDefault<bool>("syncVelocityAtBoundary",true)),
    checkDebyeLength_(readBool(particleProperties_.subDict("SolverSettings").lookup("checkDebyeLength"))),
    checkPlasmaFrequency_(readBool(particleProperties_.subDict("SolverSettings").lookup("checkPlasmaFrequency"))),
    warnCellTrajectory_(readScalar(particleProperties_.subDict("SolverSettings").lookup("warnCellTrajectory"))),
    isInitializing_(false),
    picInitialiseDict_(dictionary::null),
    constProps_(),
    rndGen_(label(149382906) + 7183*Pstream::myProcNo()),
    T_(
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    binaryCollisionModel_(),
    wallReflectionModel_(),
    maxwellSolver_(),
    particlePusher_(),
    coulombCollisionModel_(),
    electronNeutralCollisionModel_(),
    ionNeutralCollisionModel_(),
    backgroundGas_(),
    particleMerging_(),
    boundaryEvents_(*this),
    boundaryModels_(*this),
    parcelDiagnostics_(*this)
{
    //This function calculates an average cell length scale used to warnings related to the cell size being bigger than the Debye length.
    calculateCellLengthScales();

    //Read the space charge density provided in constant/picProperties
    readSpaceChargeDensity();

    //Setup the size if the cell occupancy field
    for(label i=0;i<mesh_.nCells();i++)
    {
        sortedCellOccupancy_[i].setSize(typeIdList_.size());
    }

    //Build constant parcel properties
    buildConstProps();

    //Read lagrangian fields from the time directory if needed
    if (readFields)
    {
        ParcelType::readFields(*this);
    }

    //Setup parcel properties(charge) when the typeId is known(was read)
    setupParticleProperties();

    //Setup all submodels at this point constProps are known
    setModels();

    //Build the cell occupancy
    buildCellOccupancy();

    //Inital update of the electric field
    this->eField_().update();
}

/*
Foam::PICCloud<ParcelType>::PICCloud

The constructor: Setup all members of the class.
This one is called by picInitialise.
*/
template<class ParcelType>
Foam::PICCloud<ParcelType>::PICCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& picInitialiseDict
)
    :
    Cloud<ParcelType>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.subDict("SolverSettings").lookup("typeIdList")),
    nParticleEqWeight_(0.0),
    chargedSpecies_(),
    neutralSpecies_(),
    ionSpecies_(),
    neutralTypeList_(),
    ionTypeList_(),
    electronTypeId_(-1),
    fieldCalculation_(typeIdList_.size(),-1),
    sortedCellOccupancy_(mesh_.nCells()),
    q_
    (
        IOobject
        (
            this->name() + "q_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(1, 0, -3, 0, 0), 0)
    ),
    fD_
    (
        IOobject
        (
            this->name() + "fD_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            Zero
        )
    ),
    rhoNSpecies_(),
    rhoN_
    (
        IOobject
        (
            this->name() + "rhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), vSmall)
    ),
    rhoMSpecies_(),
    rhoM_
    (
        IOobject
        (
            this->name() + "rhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), vSmall)
    ),
    picRhoN_
    (
        IOobject
        (
            this->name() + "picRhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    linearKESpecies_(),
    linearKE_
    (
        IOobject
        (
            this->name() + "linearKE_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    momentumSpecies_(),
    momentum_
    (
        IOobject
        (
            this->name() + "momentum_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -1, 0, 0),
            Zero
        )
    ),
    rhoChargeSpecies_(),
    chargeDensity_(),
    N_(),
    eField_(),
    phiE_
        (
            IOobject
            (
                "phiE",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        ),
    jSpecies_(),
    j_
    (
        IOobject
        (
            "j",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero",  dimensionSet(0, -2, 0, 0, 0,1,0), Zero),
        calculatedFvPatchScalarField::typeName
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    jE_
    (
        IOobject
        (
            "jE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -3, 0, 0,0,0), 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    spaceChargeDensity_
    (
        IOobject
        (
            "spaceChargeDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0,1,0), Zero)
    ),
    cellLengthScale_(mesh_.nCells(),0.0),
    avgcellLengthScale_(0.0),
    printVelocityWarning_(true),
    syncVelocityAtBoundary_(particleProperties_.subDict("SolverSettings").lookupOrDefault<bool>("syncVelocityAtBoundary",true)),
    checkDebyeLength_(readBool(particleProperties_.subDict("SolverSettings").lookup("checkDebyeLength"))),
    checkPlasmaFrequency_(readBool(particleProperties_.subDict("SolverSettings").lookup("checkPlasmaFrequency"))),
    warnCellTrajectory_(readScalar(particleProperties_.subDict("SolverSettings").lookup("warnCellTrajectory"))),
    isInitializing_(true),
    picInitialiseDict_(picInitialiseDict),
    constProps_(),
    rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 0, 1, 0), 0.0)
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -1, 0, 0),
            Zero
        )
    ),
    binaryCollisionModel_(),
    wallReflectionModel_(),
    maxwellSolver_(),
    particlePusher_(),
    coulombCollisionModel_(),
    electronNeutralCollisionModel_(),
    ionNeutralCollisionModel_(),
    backgroundGas_(),
    particleMerging_(),
    boundaryEvents_(*this),
    boundaryModels_(*this),
    parcelDiagnostics_(*this)
{

    //This function calculates an average cell length scale used for warnings related to the cell size being bigger than the Debye length.
    calculateCellLengthScales();

    //Read the space charge density provided in constant/picProperties
    readSpaceChargeDensity();

    //Remove already existing particles?
    bool clearParts = readBool(picInitialiseDict.lookup("clearParticles"));
    if(clearParts) {
        clear();
    }
    else
        ParcelType::readFields(*this);

    //Setup the size if the cell occupancy field
    for(label i=0;i<mesh_.nCells();i++)
    {
        sortedCellOccupancy_[i].setSize(typeIdList_.size());
    }

    //Build constant parcel properties
    buildConstProps();

    //Setup parcel properties(charge) when the typeId is known(was read)
    setupParticleProperties();

    //Setup all submodels at this point constProps are known
    setModels();

    //Inital update of the electric field
    this->eField_().update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PICCloud<ParcelType>::~PICCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
Foam::PICCloud<ParcelType>::calculateCellLengthScales

Calculate the cell length scale. Used to warn the used if the cell is bigger than the Debye length.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::calculateCellLengthScales()
{
    //In 3D simulation this is simply the cube root of the cells volume
    cellLengthScale_ = mag(cbrt(mesh_.V()));

    //In 1D,2D ignore invalid geometric directions
    if(mesh_.nSolutionD() < 3)//There is probably a better way calculating this...
    {
        const volVectorField& C = mesh_.C();
        const surfaceVectorField& Cf = mesh_.Cf();

        vectorField cellLength (C.size(),Zero);

        //Get valid coords
        vector validCoords(1.0,1.0,1.0);
        meshTools::constrainDirection(mesh_, mesh_.solutionD(), validCoords);

        //For all cells...
        forAll(mesh_.cells(),celli)
        {
            cell c = mesh_.cells()[celli];
            //and for all the faces of the cell...
            forAll(c,fi)
            {

                label facei = c[fi];

                if(facei >= Cf.size())//only internal faces
                    continue;

                //add the distance to the cells face
                cellLength[celli] += vector(
                                mag(Cf[facei].x()-C[celli].x()),
                                mag(Cf[facei].y()-C[celli].y()),
                                mag(Cf[facei].z()-C[celli].z())
                            );
            }

        }

        //Add the distance to the boundary faces, which are saved in a different list.
        forAll(mesh_.boundary(),patchi)
        {
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            const vectorField& pCf = Cf.boundaryField()[patchi];
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                cellLength[own] += vector(
                                mag(pCf[pFacei].x()-C[own].x()),
                                mag(pCf[pFacei].y()-C[own].y()),
                                mag(pCf[pFacei].z()-C[own].z())
                            );
            }
        }

        //Set the cell length in invalid directions to zero, for the valid directions calculate the absolut value directly (1D) or calculate the length scale the via the area (2D)
        forAll(cellLength,i) {
                cellLength[i].x() *= validCoords.x();
                cellLength[i].y() *= validCoords.y();
                cellLength[i].z() *= validCoords.z();

                if(mesh_.nSolutionD() == 1) {
                    cellLengthScale_[i] = mag(cellLength[i]);
                }
                else {
                    cellLength[i] += (-validCoords+vector(1.0,1.0,1.0));//Set invalid dirs to one so we do not mulitply by zero
                    scalar area = cellLength[i].x()*cellLength[i].y()*cellLength[i].z();
                    cellLengthScale_[i] = sqrt(area);
                }
        }
    }
    //Calculate an average length scale and communicate for parallel runs
    avgcellLengthScale_ = average(cellLengthScale_);
    if(Pstream::parRun())
    {
        reduce(avgcellLengthScale_,sumOp<scalar>());
        avgcellLengthScale_ /= Pstream::nProcs();
    }
}

/*
Foam::PICCloud<ParcelType>::readSpaceChargeDensity()

Read the space charge and calculate the space charge density field.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::readSpaceChargeDensity()
{
    //Get the number of cells in the mesh
    label cellCount = mesh_.cells().size();
    reduce(cellCount,sumOp<label>());

    //Read value in coulombs
    scalar spCh = readScalar(particleProperties_.subDict("SolverSettings").lookup("spaceCharge"));
    Info << nl << "Read space charge value of " << spCh << " As" << endl;

    //The space charge is distributed uniformly to all cells...
    spCh/=scalar(cellCount);

    forAll(mesh_.cells(), cellI)
    {
        spaceChargeDensity_[cellI] = spCh;
    }

    scalarField& spaceChargeDensity = spaceChargeDensity_.primitiveFieldRef();
    spaceChargeDensity/=mesh_.cellVolumes();
}

/*
Foam::PICCloud<ParcelType>::evolve()

Evolve the PIC cloud. This one is called by picFoam every time step.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::evolve()
{
    //Tracking data. Is to keep track of properties during the parcels movement.
    typename ParcelType::trackingData td(*this);

    //Reset fields.
    resetFields();

    //Call the boundary model to inject new parcels into the domian.
    boundaryModels().injection();

    //If the velocity of one parcel is too high, warn once on each processor. After this, set this variable too false, so we do not spam stdout.
    printVelocityWarning_ = true;
    //Moves the parcels. Calls into the Cloud base class which then calls back to Foam::PICParcel<ParcelType>::move.
    Cloud<ParcelType>::move(*this, td, mesh_.time().deltaTValue());

    //Call the boundary event model after all parcels have been moved.
    boundaryEventModels().postMove();

    //Build the cell occupancy.
    buildCellOccupancy();

    //Update the background gas.
    backgroundGas().update();

    //Call all collision algorithm.
    collisions();

    //Calculate the fields.
    calculateFields();

    //Update boundary models before the new electric field is calculated (used by some circuit models).
    boundaryModels().preUpdate_Boundary();

    //Solve the electric field.
    maxwellSolver().solveFields();

    //Update boundary models after the electric field has been updated.
    boundaryModels().postUpdate_Boundary();

    //Update joule heat
    jE_ = j_ & eField_().field();

    //Merge parcels
    particleMerging().checkAndMerge();
}


/*
Foam::PICCloud<ParcelType>::info()

Calculated and print diagnostics.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::info()
{
    //Print boundary event info.
    boundaryEventModels().info();

    //Get number the of parcels.
    label nPICParticles = this->size();
    reduce(nPICParticles, sumOp<label>());

    //Display global Info
    Info << nl << "[Diagnostics]" << nl
              << "Cloud name: " << this->name() << nl
              << "    Number of pic particles        = "
              << nPICParticles
              << nl << "----" << endl;//flush

    //If we have parcels print diagnostics.
    if(nPICParticles)
        Diagnostics().info();

    //During the simulation, if required warn if the cell is larger than the average debye length.
    if(checkDebyeLength_ && !isInitializing())
    {
        //Calculate the debye length.
        tmp<volScalarField::Internal> debye = debyeField();

        //Get the average value.
        scalar avgDebyeL = average(debye.ref().field());
        if(Pstream::parRun())
        {
            //Average over all processor in a parallel run.
            reduce(avgDebyeL,sumOp<scalar>());
            avgDebyeL /= Pstream::nProcs();
        }

        Info << "    Average cell size: " << avgcellLengthScale_ << " m^-3" << nl
             << "    Average debye length: " << avgDebyeL << " m^-3" << endl;
        if(avgcellLengthScale_ > avgDebyeL)
        {
            Info << "    => Average cell size is larger than the average debye length!" << endl;
        }
        //Write the debye field if it's time for this.
        if (this->db().time().writeTime())
            debye.ref().write();
    }
    
    //During the simulation, if required warn if the time step is to big for the current plasma frequency.
    if(checkPlasmaFrequency_ && !isInitializing())
    {
        //Calculate the plasma frequency.
        tmp<volScalarField::Internal> plasmaFreq = plasmaFreqField();

        //Get the average...
        scalar avgPlasmaFreq = average(plasmaFreq.ref().field());
        if(Pstream::parRun())
        {
            //Average over all processor in parallel run.
            reduce(avgPlasmaFreq,sumOp<scalar>());
            avgPlasmaFreq /= Pstream::nProcs();
        }
        Info << "    Average plasma frequency: " << avgPlasmaFreq  << " Hz" << endl;
        if(avgPlasmaFreq*mesh_.time().deltaTValue() > 2.0)
        {
            Info << "    => Average plasma frequency times deltaT is larger than 2." << endl;
        }
    }
}

/*
Foam::PICCloud<ParcelType>::equipartitionLinearVelocity()

Sample a velocity according to the Maxwell-Boltzmann distribution.
*/
template<class ParcelType>
Foam::vector Foam::PICCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.sampleNormal<vector>();
}

/*
Foam::PICCloud<ParcelType>::dumpParticle()

Dump parcels to a file.
*/
template<class ParcelType>
void Foam::PICCloud<ParcelType>::dumpParticle() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcels_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(typename PICCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        pObj<< "x " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << " u " << p.U().x()
            << " "  << p.U().y()
            << " "  << p.U().y()
            << " t "  << p.typeId()
            << " n " << p.nParticle()
            << " c " << p.chargeModifier()
            << nl;
    }

    pObj.flush();
}


template<class ParcelType>
void Foam::PICCloud<ParcelType>::autoMap(const mapPolyMesh& mapper)
{
    //Usage of dynamic meshes is not supported
    FatalErrorInFunction << "Dynamic meshes are not supported" << abort(FatalError);
    //Cloud<ParcelType>::autoMap(mapper);
}

/*
Foam::PICCloud<ParcelType>::plasmaFreqField()
Calculate the plasma frequency
*/
template<class ParcelType>
tmp<volScalarField::Internal> Foam::PICCloud<ParcelType>::plasmaFreqField() const
{
   tmp<volScalarField::Internal> tplasmaFrequency(
       volScalarField::Internal::New
       (
            "plasmaFrequency",
             mesh_,
             dimensionedScalar("zero",  dimless/dimTime, 0.0)
       ));

   Field<scalar>& plasmaFrequency = tplasmaFrequency.ref().field();
   const typename ParcelType::constantProperties& cP = constProps(electronTypeId_);
   const List<List<DynamicList<ParcelType*>>>& sortedCellOccupancy(this->sortedCellOccupancy());

   //Go through all cells...
   forAll(mesh_.cells(),celli)
   {
       //Does not include new ionized parcels if called after collision...
       forAll(sortedCellOccupancy[celli][electronTypeId_],parti)//Only electrons...
       {
           ParcelType* p = sortedCellOccupancy[celli][electronTypeId_][parti];
           if(p == nullptr)//sortedCellOccupancy wasn't re-evaluated some parcels are removed e.g. by ionization code.
               continue;
           plasmaFrequency[celli] += p->nParticle();//Add the parcel weight.
       }
       plasmaFrequency[celli] /= mesh_.cellVolumes()[celli];//Divide by the volume => Get the number density.
   }
   plasmaFrequency = sqrt(plasmaFrequency*cP.charge()*cP.charge()/constant::electromagnetic::epsilon0.value()/cP.mass());//Calculate the plasma frequency
   return tplasmaFrequency;
}

/*
Foam::PICCloud<ParcelType>::debyeField()
Calculate the field containing the Debye length for each cell.
*/
template<class ParcelType>
tmp<volScalarField::Internal> Foam::PICCloud<ParcelType>::debyeField() const
{
    tmp<volScalarField::Internal> tdebyeField(
      volScalarField::Internal::New
      (
          "debyeLength",
           mesh_,
           dimensionedScalar("zero",  dimLength, 0.0)
      ));

      volScalarField::Internal& debyeField = tdebyeField.ref();

        const List<List<DynamicList<ParcelType*>>>& sortedCellOccupancy(this->sortedCellOccupancy());
        //For all cells...
        forAll(mesh_.cells(),celli)
        {
            scalar debye(0.0);
            //For all charged species...
            forAll(chargedSpecies(),si)
            {
                scalar  T(0.0), n(0.0);
                label speci = chargedSpecies()[si];
                scalar mass = constProps(speci).mass();

                scalar charge = 0.0;
                scalar vSqr = 0.0;
                vector v = Zero;

                //Does not include new ionized parcels if called after collision...
                forAll(sortedCellOccupancy[celli][speci],parti)
                {
                    ParcelType* p = sortedCellOccupancy[celli][speci][parti];
                    if(p == nullptr)//ignore removed parcels
                        continue;

                    vSqr += (p->U() & p->U())*p->nParticle();
                    charge += p->charge()*p->nParticle();
                    v += p->U()*p->nParticle();
                    n += p->nParticle();
                }
                if(n <= 0.0) {
                    continue;
                }
                //Average velocity (drift velocity).
                v /= n;
                //Average charge.
                charge /= n;
                //Calculate the average squared velocity.
                vSqr /= n;
                //Number density of the species.
                n /= mesh().cellVolumes()[celli];

                //Calculate the temperature based on the Maxwell-Boltzmann distribution.
                T = mass/(3.0*constant::physicoChemical::k.value())*(vSqr-(v&v));

                if(T <= 0.0) {
                    continue;
                }

                //Calculate the Debye length.
                debye += n*charge*charge/constant::electromagnetic::epsilon0.value() * 1.0/(T*constant::physicoChemical::k.value());
            }
            if(debye > 0.0)//If we have a valid value update the cell value of the field.
                debyeField[celli] = ::sqrt(1/debye);
        }
        return tdebyeField;
}

// ************************************************************************* //
