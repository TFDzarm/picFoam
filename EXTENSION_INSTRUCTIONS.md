# Instructions for extending the code
***picFoam*** is build in a modular way, the main application is included in the *solver/picFoam* directory.
The preprocessing tool ***picInitalise*** is included in the *preProcess/picInitialise* directory.
Both build upon the PIC-MCC libary included in the *PICMCCFVM* directory.
 
## Structur of the PICMCCFVM libary

The librarie's directory contains three additional directories:

1. ***clouds***: The code in this directory definies the class PICCloud (see: *clouds/Templates/PICCloud/PICCloud.H*) which at its core is dereived from a double linked list
containing the parcels created at runtime. Class members of PICCloud include smart pointers to the submodels choosen at runtime, fields and all functions needed for the PIC algorithm.
The most important function is PICCloud::evolve called by the solver ***picFoam*** in its main loop.
2. ***parcels***: The code in this directory defines the class PICParcel (see: *parcels/Templates/PICParcel/PICParcel.H*). It represents a number of real particles (electrons, atoms, molecules).
Functions of this class definied in the source file PICParcel.C include calls to the ParticlePusher submodule to update the parcel's velocity. Calls to the Particle base class for moving the parcel inside the mesh
and calls to the boundary models to handle the interaction with them.
3. ***submodels***: This directory contains subdirectories for all submodels defined in ***picFoam***.

These include:

* **BackgroundGasModel**: Neutral background gas, defined via a density and temperature fields and interacting with particles in the collision algorithm.
* **BoundaryEvent**: Models that are called when a particle hits a boundary (mainly for diagnostic purpose).
* **BoundaryModel**: Boundary models which are responsible for injecting new particles and defining the boundary interaction (e.g. circuits). 
* **ChargeDistribution**: Models which distribute the particle's charge to the charge density field *rhoCharge*.
* ***CollisionModels***
	* **BinaryCollisionModel**: This model handles binary collisions between neutral species.
		* **TotalCrossSectionModel**: Models which define the total cross section for species undergoing binary collision.
	* **CoulombCollisionModel**: This model handles coulomb collisions between charged particles.
		* **PairingAlgorithm**: Models which are responsible for pairing charged particle to perform a single binary collision with a large scattering angle.
		* **IonizationModel**: Models that describe the ionization occurring in the collision of charged particles.
	* **ElectronNeutralCollisionModel**: Collisions between electrons and neutral species.
	* **IonNeutralCollisionModel**: Collisions between one ion and one neutral species.
	* **CrossSectionModel**: Models defining the cross section for different interactions between electrons and neutrals as well as ion and neutrals.
	* **WeightCorrectionModel**: Models for correcting the post-collision velocities in the case of a collision between parcels of unequal weight.
* **Diagnostics**: Diagnostic calculations.
* **FieldWeighting**: Models for weighing the electric field to the particle's position.
* **MaxwellSolver**: Models for the calculation of the electric field.
* **ParticleMerging**: Particle merging to reduce the computational cost of the simulation.
* **ParticlePusher**: Algorithms to update the parcel's velocity in electromagnetic fields.
* **WallReflectionModel**: Models for interactions with solid walls.

## How to add a new submodel?

The addition of a new model includes two steps:

* To add a new submodel to one of the submodel base classes, listed above, one has to inherit from the base class interface for the desired model. These are defined in the corresponding subdiretory e.g. *submodels/ParticlePusher/ParticlePusher*. The easiest way of doing this is to copy an existing model and to rename it.
Hereby the user has to make sure to define a unique name with the TypeName macro inside the class declaration (see e.g. line 62 in *submodels/ParticlePusher/BorisPusher/BorisPusher.H*). The name definined with this macro is saved in an *RunTimeSelectionTable* from where the model is choosen via the picProperties file at runtime.

* Once the code for the new model is written the new model has to be exposed. For this, the user has to add a line to corresponding submodel file in the directory **parcels/derived/picParcel/...**.

For example to add an imaginary ParticlePusher model called MyPusher the user adds the class name via the appropriate macro:

		...
        //Add this line:
		#include "MyPusher.H"
		...
		namespace Foam
		{
			typedef PICCloud<picParcel> CloudType;

			makeParticlePusher(PICCloud<picParcel>)

			makeParticlePusherType(BorisPusher, CloudType)
    		makeParticlePusherType(BorisNRPusher, CloudType)
    		makeParticlePusherType(VayPusher, CloudType)
    		makeParticlePusherType(HigueraCaryPusher, CloudType)
			//Add this line:
			makeParticlePusherType(MyPusher, CloudType)
		}

Since all submodels are templated C++ classes, the definition of the model's functions (the code in the source file e.g. MyPusher.C) need to be included with the class declaration. This is because templates are instantiated at compile-time not link-time! For this reason, at the end of the model's header file each model includes the source file:

		#ifdef NoRepository
			#include "MyPusher.C"
		#endif
		// Such as at line 95 in submodels/ParticlePusher/BorisPusher/BorisPusher.H

This in a way makes each model header-only and the compiler does not need to be informed about the new files. If one adds non-templated code to the library the path to the source file has to be included in ***Make/files***. Finally, after adding the new code, the user has to run the command *wmake* in the root directory of the library to recompile it, the updated library is automaticlly added to the user directory of OpenFOAM.
