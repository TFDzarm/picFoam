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

Class
    Foam::PICCloud

Description
    Templated base class for pic cloud

SourceFiles
    PICCloudI.H
    PICCloud.C

\*---------------------------------------------------------------------------*/

#ifndef PICCloud_H
#define PICCloud_H

#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "barycentric.H"

#include "pointMesh.H"
#include "pointBoundaryMesh.H"
#include "pointFields.H"

#include "FieldWeigthing.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class BinaryCollisionModel;

template<class CloudType>
class WallReflectionModel;

template<class CloudType>
class BoundaryModel;

//PIC Submodels
template<class CloudType>
class MaxwellSolver;

template<class CloudType>
class ParticlePusher;

template<class CloudType>
class CoulombCollisionModel;

template<class CloudType>
class BoundaryEventModelList;

template<class CloudType>
class BoundaryModelList;

template<class CloudType>
class DiagnosticsList;

template<class CloudType>
class ElectronNeutralCollisionModel;

template<class CloudType>
class IonNeutralCollisionModel;

template<class CloudType>
class ChargeDistribution;

template<class CloudType>
class BackgroundGasModel;

template<class CloudType>
class ParticleMerging;

//---------------------------//

/*---------------------------------------------------------------------------*\
                         Class PICCloud Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PICCloud);

template<class ParcelType>
class PICCloud
:
    public Cloud<ParcelType>,
    public PICCloudName
{
    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudName_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of instances of species types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (Ar Ne+ e) Ar has typeId label = 0, Ne+ = 1, e = 2.
        List<word> typeIdList_;

        //- Number of real atoms represented by a parcel (on equal weight / else == 0.0)
        scalar nParticleEqWeight_;

        //List holding the typeId of all species holding a charge
        List<label> chargedSpecies_;


        //Setup by SpeciesRelations:
        // -List holding the typeId of all species that are neutral
        List<label> neutralSpecies_;

        // -List holding the typeId of all species that are ions
        List<label> ionSpecies_;

        // -List holding the typeId of the neutral/atomic species, -1 if not present
        List<label> neutralTypeList_;

        // -List holding the typeId of the ion species, -1 if not present
        List<label> ionTypeList_;

        // -Typeid of the electron species
        label electronTypeId_;


        //index of species into field lists; -1 no calculation
        List<label> fieldCalculation_;

        //A data structure holding which particles are in which cell (sorted by typeId)
        //- 1st list cells
        //- 2nd list typeIds
        //- 3rd list parcel pointer
        List<List<DynamicList<ParcelType*>>> sortedCellOccupancy_;

        //- Heat flux at surface field
        volScalarField q_;

        //- Force density at surface field
        volVectorField fD_;

        //- Number density fields
        PtrList<volScalarField> rhoNSpecies_;
        volScalarField rhoN_;

        //- Mass density fields
        PtrList<volScalarField> rhoMSpecies_;
        volScalarField rhoM_;

        //- PIC particle density field
        volScalarField picRhoN_;

        //- Linear kinetic energy density fields
        PtrList<volScalarField> linearKESpecies_;
        volScalarField linearKE_;

        //- Momentum density fields
        PtrList<volVectorField> momentumSpecies_;
        volVectorField momentum_;

        // Charge density fields
        PtrList<volScalarField> rhoChargeSpecies_;
        autoPtr<ChargeDistribution<PICCloud<ParcelType>>> chargeDensity_;

        //Density fields
        PtrList<volScalarField> N_;

        //Electric field model
        autoPtr<FieldWeigthing> eField_;

        //Electric potential field
        volScalarField phiE_;

        //Current density fields
        PtrList<volVectorField> jSpecies_;
        volVectorField j_;

        //Magnetic field
        volVectorField B_;

        //Joule heat
        volScalarField jE_;

        //Space charge density field
        volScalarField spaceChargeDensity_;

        //Cell length scales
        scalarField cellLengthScale_;
        scalar avgcellLengthScale_;

        //Used to suppress the printing of multiple warings per thread
        bool printVelocityWarning_;

        //Synchronise leap frog velocity on boundary collision (default=true)
        bool syncVelocityAtBoundary_;

        //Check cell size against debye length every timestep
        bool checkDebyeLength_;

        //Check dt against plasma frequency every timestep
        bool checkPlasmaFrequency_;

        //Warn if particle moved more then this portion of the cell length scale
        scalar warnCellTrajectory_;

        //Do we run the picFoam or picInitialize?
        bool isInitializing_;

        //picInitialiseDict, null if not initialising
        dictionary picInitialiseDict_;

        //- Parcel constant properties - one for each type
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        // boundary value fields

            //- Boundary temperature
            volScalarField T_;

            //- Boundary velocity
            volVectorField U_;


        // References to the cloud sub-models

            //- Binary (Neutral + Ion and Neutral) collision model
            autoPtr<BinaryCollisionModel<PICCloud<ParcelType>>>
                binaryCollisionModel_;

            //- Wall interaction model
            autoPtr<WallReflectionModel<PICCloud<ParcelType>>>
                wallReflectionModel_;

            //- Maxwell field solver
            autoPtr<MaxwellSolver<PICCloud<ParcelType>>>
                maxwellSolver_;

            //- Particle pusher/Velocity integrator
            autoPtr<ParticlePusher<PICCloud<ParcelType>>>
                particlePusher_;

            //- Coulomb collision model
            autoPtr<CoulombCollisionModel<PICCloud<ParcelType>>>
                coulombCollisionModel_;

            //- Electron - Neutral collision model
            autoPtr<ElectronNeutralCollisionModel<PICCloud<ParcelType>>>
                electronNeutralCollisionModel_;

            //- Ion - Neutral collision model
            autoPtr<IonNeutralCollisionModel<PICCloud<ParcelType>>>
                ionNeutralCollisionModel_;

            //- Background gas model
            autoPtr<BackgroundGasModel<PICCloud<ParcelType>>>
                backgroundGas_;

            //- Particle merging model
            autoPtr<ParticleMerging<PICCloud<ParcelType>>>
                particleMerging_;

    // Private Member Functions

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Setup non-const particle properties
        void setupParticleProperties();

        //- Calculate collisions between molecules
        void collisions();

        //- Reset the data accumulation field values to zero
        void resetFields();

        //- Calculate cell length scales
        void calculateCellLengthScales();

        //- Read space charge from dictionary and calculate the density
        void readSpaceChargeDensity();
private:
        //- Disallow default bitwise copy construct
        PICCloud(const PICCloud&);

        //- Disallow default bitwise assignment
        void operator=(const PICCloud&);


public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //  file
        PICCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        PICCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& picInitialiseDict
        );


    //- Destructor
    virtual ~PICCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudName() const;

                //- Return references to the mesh
                inline const fvMesh& mesh() const;

                //- Calculate debye length field
                tmp<volScalarField::Internal> debyeField() const;

                //- Calculate plasma frequency field
                tmp<volScalarField::Internal> plasmaFreqField() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return the idList
                inline const List<word>& typeIdList() const;

                //- Return species fieldCalculation fields
                inline const List<label>& fieldCalculation() const;

                //- Return the idList of charged species
                inline const List<label>& chargedSpecies() const;

                //- Return the idList of neutral species
                inline const List<label>& neutralSpecies() const;

                //- Return the idList of ion species
                inline const List<label>& ionSpecies() const;

                //- Return the number of real particles represented by one
                //  parcel (on equal weight / else == 0.0)
                inline scalar nParticleEqWeight() const;

                //- Return the sorted cell occupancy addressing
                inline const List<List<DynamicList<ParcelType*>>>&
                    sortedCellOccupancy() const;

                inline List<List<DynamicList<ParcelType*>>>&
                    sortedCellOccupancy();

                //- Return all of the constant properties
                inline const List<typename ParcelType::constantProperties>&
                    constProps() const;

                //- Return the constant properties of the given typeId
                inline const typename ParcelType::constantProperties&
                    constProps(label typeId) const;

                //- Return references to the random object
                inline Random& rndGen();

            // References to the boundary fields for surface data collection

                //- Return non-const heat flux boundary field reference
                inline volScalarField::Boundary& qBF();

               //- Return non-const force density at boundary field reference
               inline volVectorField::Boundary& fDBF();

                //- Return non-const number density boundary field reference
                inline volScalarField::Boundary& rhoNSpeciesBF(label fieldIndex);
                inline volScalarField::Boundary& rhoNBF();

                //- Return non-const mass density boundary field reference
                inline volScalarField::Boundary& rhoMSpeciesBF(label fieldIndex);
                inline volScalarField::Boundary& rhoMBF();

                //- Return non-const linear kinetic energy density boundary
                //  field reference
                inline volScalarField::Boundary& linearKESpeciesBF(label fieldIndex);
                inline volScalarField::Boundary& linearKEBF();

                //- Return non-const momentum density boundary field reference
                inline volVectorField::Boundary& momentumSpeciesBF(label fieldIndex);
                inline volVectorField::Boundary& momentumBF();

                inline volVectorField::Boundary& jBF();

            // References to the macroscopic fields

                //- Return heat flux at surface field
                inline const volScalarField& q() const;

                //- Return force density at surface field
                inline const volVectorField& fD() const;

                //- Return macroscopic temperature - boundary values only
                inline const volScalarField& T() const;

                //- Return macroscopic velocity - boundary values only
                inline const volVectorField& U() const;

                //- Return the real particle number density field
                inline const volScalarField& rhoN() const;

                //- Return the particle mass density field
                inline const volScalarField& rhoM() const;

                //- Return the field of number of PIC particles
                inline const volScalarField& picRhoN() const;

                //- Return the total linear kinetic energy (translational and
                // thermal density field
                inline const volScalarField& linearKE() const;

                //- Return the momentum density field
                inline const volVectorField& momentum() const;

               // Non-const access
               //- Return non-const reference to the charge density field
               inline volScalarField& rhoCharge();

               //- Return non-const reference to the electric field
               inline volVectorField& electricField();

               //- Return non-const reference to the electric potentialField
               inline volScalarField& elpotentialField();

               //- Return non-const reference to the magentic field
               inline volVectorField& magneticField();

                //Const access
               //- Return const reference to the charge density field
               inline const volScalarField& rhoCharge() const;

               //- Return const reference to the electric field
               inline const volVectorField& electricField() const;

               //- Return field weighting model
               inline FieldWeigthing& eFieldWeighting();
               inline const FieldWeigthing& eFieldWeighting() const;

               //- Return non-const reference to the electric potentialField
               inline const volScalarField& elpotentialField() const;

               //- Return const reference to the magentic field
               inline const volVectorField& magneticField() const;

               //- Return const reference to the space charge density
               inline const volScalarField& spaceChargeDensity() const;

               //- Do we synchronise the leap frog velocity on boundary collision?
               inline const bool& syncVelocityAtBoundary() const;

               //- Return non-const access to the joule heat
               inline volScalarField& jouleHeat();

               //- Return const access to the joule heat
               inline const volScalarField& jouleHeat() const;

               //- Return the cell length scales
               inline const scalarField& cellLengthScale() const;
               inline const scalar& avgCellLengthScale() const;

               //- Do we still print a waring in this thead?
               inline bool& printVelocityWarning();
               inline const bool& printVelocityWarning() const;

               //- Return the electron species type id
               inline const label& electronTypeId() const;

               //- Return the neutral species for typeId i (-1 if it does not exist)
               inline const label& neutralTypeId(label i) const;

               //- Return ion typeId list
               inline const List<label>& ionTypeIdList() const;

               //- Return the ion species for typeId i (-1 if it does not exist)
               inline const label& ionTypeId(label i) const;

               //- Return the fraction of the cell length scale before we print a warning
               inline const scalar& warnCellTrajectory() const;

            // Kinetic theory helper functions

                //- Generate a random velocity sampled from the Maxwellian speed
                // distribution
                vector equipartitionLinearVelocity
                (
                    scalar temperature,
                    scalar mass
                );

                // From the Maxwellian distribution:
                //- Average particle speed
                inline scalar maxwellianAverageSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianAverageSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- RootMeanSquare particle speed
                inline scalar maxwellianRMSSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianRMSSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- Most probable speed
                inline scalar maxwellianMostProbableSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianMostProbableSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                inline scalar maxwellianComponent(
                        scalar temperature,
                        scalar mass
                ) const;

                inline scalar lorentzFactor(
                        scalar velocity
                ) const;

                inline scalar lorentzFactor(
                        scalar momentum,
                        scalar mass
                ) const;

            // Sub-models

                //- Return reference to binary elastic collision model
                inline const BinaryCollisionModel<PICCloud<ParcelType>>&
                    binaryCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline BinaryCollisionModel<PICCloud<ParcelType>>&
                    binaryCollision();

                //- Return reference to wall interaction model
                inline const WallReflectionModel<PICCloud<ParcelType>>&
                    wallReflection() const;

                //- Return non-const reference to wall interaction model
                inline WallReflectionModel<PICCloud<ParcelType>>&
                    wallReflection();

                //- Return reference to MaxwellSolver
                inline const MaxwellSolver<PICCloud<ParcelType>>&
                    maxwellSolver() const;

                //- Return non-const reference to MaxwellSolver
                inline MaxwellSolver<PICCloud<ParcelType>>&
                    maxwellSolver();

                //- Return reference to ParticlePusher
                inline const ParticlePusher<PICCloud<ParcelType>>&
                    particlePusher() const;

                //- Return non-const reference to ParticlePusher
                inline ParticlePusher<PICCloud<ParcelType>>&
                    particlePusher();

                //- Return reference to ParticleMerging
                inline const ParticleMerging<PICCloud<ParcelType>>&
                    particleMerging() const;

                //- Return non-const reference to ParticleMerging
                inline ParticleMerging<PICCloud<ParcelType>>&
                    particleMerging();

                //- Return reference to CoulombCollisionModel
                inline const CoulombCollisionModel<PICCloud<ParcelType>>&
                    coulombCollision() const;

                //- Return non-const reference to CoulombCollisionModel
                inline CoulombCollisionModel<PICCloud<ParcelType>>&
                    coulombCollision();

                //- Return reference to binary elastic collision model
                inline const ElectronNeutralCollisionModel<PICCloud<ParcelType>>&
                    electronNeutralCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline ElectronNeutralCollisionModel<PICCloud<ParcelType>>&
                    electronNeutralCollision();

                //- Return reference to the ion-neutral collision model
                inline const IonNeutralCollisionModel<PICCloud<ParcelType>>&
                ionNeutralCollision() const;

                //- Return non-const reference to the ion-neutral collision model
                inline IonNeutralCollisionModel<PICCloud<ParcelType>>&
                ionNeutralCollision();

                //- Return non-const reference to BackgroundGasModel
                inline BackgroundGasModel<PICCloud<ParcelType>>&
                    backgroundGas();

                //- Return reference to BackgroundGasModel
                inline const BackgroundGasModel<PICCloud<ParcelType>>&
                    backgroundGas() const;

private:
//Define after parcelType
                //- Boundary event model list
                BoundaryEventModelList<PICCloud<ParcelType>> boundaryEvents_;

                //- Boundary model list
                BoundaryModelList<PICCloud<ParcelType>> boundaryModels_;

                //- Diagnostics list
                DiagnosticsList<PICCloud<ParcelType>> parcelDiagnostics_;
public:
                //- Return reference to BoundaryEventModelList
                inline const BoundaryEventModelList<PICCloud<ParcelType>>&
                    boundaryEventModels() const;

                //- Return non-const reference to BoundaryEventModelList
                inline BoundaryEventModelList<PICCloud<ParcelType>>&
                    boundaryEventModels();

                //- Return reference to BoundaryModelList
                inline const BoundaryModelList<PICCloud<ParcelType>>&
                    boundaryModels() const;

                //- Return non-const reference to BoundaryModelList
                inline BoundaryModelList<PICCloud<ParcelType>>&
                    boundaryModels();

                //- Return reference to DiagnosticsList
                inline const DiagnosticsList<PICCloud<ParcelType>>&
                    Diagnostics() const;

                //- Return non-const reference to DiagnosticsList
                inline DiagnosticsList<PICCloud<ParcelType>>&
                    Diagnostics();

        // Check

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total internal energy in the system
            inline scalar internalEnergyOfSystem() const;

            //- Print cloud information / diagnostics
            void info();

            //- Dump particle positions to .obj file
            void dumpParticle() const;


            void setModels();

            inline const bool& isInitializing() const;

            inline dictionary& picInitialiseDict();

            inline const dictionary& picInitialiseDict() const;

        // Cloud evolution functions

            //- Add new parcel
            ParcelType* addNewParcel
            (
                const vector& position,
                const label celli,
                const vector& U,
                const label typeId
            );

            //- Record which particles are in which cell
            void buildCellOccupancy();

            //- Calculate the volume field data
            void calculateFields();

            //- Evolve the cloud (move, collide, ...)
            void evolve();

            //- Clear the Cloud
            inline void clear();


        // Mapping

            //- Remap the particles to the correct cells following mesh change
            // NOT IMPLEMENTED!
            virtual void autoMap(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PICCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PICCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
