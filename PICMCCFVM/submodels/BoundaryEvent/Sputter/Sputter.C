/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "Sputter.H"
#include "DynamicList.H"
#include "polyPatch.H"
#include "ListOps.H"
#include "constants.H"
#include "Time.H"
#include "polyMesh.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SputterEvent<CloudType>::SputterEvent
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryEvent<CloudType>(dict,cloud,typeName,associatedPatches),
    sputterSpecies_(),
    initalEnergies_(),
    sputterProbability_(),
    patchNormal_(),
    patchTriFaces_(),
    printInfo_(this->coeffDict().lookupOrDefault("printInfo",false)),
    sputterdParcels_(),
    totalSputterdParcels_(),
    runningAverage_(),
    averageStart_(cloud.mesh().time().value())
{
    label nSpecies = this->owner().typeIdList().size();

    //Get number of patches not counting processorPatches
    label nPatches = cloud.mesh().boundaryMesh().size();//FIXME: Used in multiple classes maybe calculate it globally once...
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

    sputterSpecies_.setSize(nPatches,List<List<label>>(nSpecies));
    sputterProbability_.setSize(nPatches,List<List<scalar>>(nSpecies,List<scalar>(nSpecies,0.0)));
    initalEnergies_.setSize(nPatches,List<List<scalar>>(nSpecies));

    //Read options for all patches...
    forAll(associatedPatches,i)
    {
        label patchId = associatedPatches[i];
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patchId];
        const dictionary& patchDict = this->coeffDict().subDict(patch.name());

        forAllConstIter(IDLList<entry>, patchDict, iter)
        {
            if(iter().isDict())
            {
                const dictionary& subDict = iter().dict();

                label id = findIndex(this->owner().typeIdList(),iter().keyword());//Species for which we sputter
                wordList typeIdWords = subDict.lookup("sputterTypeIds");//Species which will be sputtered/injected
                scalarList typeIdEnergies = subDict.lookup("initialEnergies");//Energie of the injected species
                initalEnergies_[patchId][id].transfer(typeIdEnergies);
                scalarList probability = subDict.lookup("probabilities");//Probabilities for the event to occur

                //Lists all need to be the same size
                if(typeIdWords.size() != initalEnergies_[patchId][id].size() && probability.size() != typeIdWords.size())
                    FatalErrorInFunction << "Number of species to be sputtered (" <<  typeIdWords.size() << ") does not equal number of inital energies (" << initalEnergies_[patchId][id].size() << ") or probabilities (" << probability.size() << ") definied" << abort(FatalError);

                DynamicList<label> sputterSpecies;
                forAll(typeIdWords,i)
                {
                    label sId = findIndex(this->owner().typeIdList(), typeIdWords[i]);
                    if(sId == -1)
                        FatalErrorInFunction << "Undefined typeId " << iter().keyword() << abort(FatalError);


                    sputterSpecies.append(sId);
                }
                sputterSpecies_[patchId][id].transfer(sputterSpecies);
                sputterProbability_[patchId][id].transfer(probability);

            }

        }
    }

    //Setup list for the injection code
    //Get the triangles that make up a faces, used for calculating an in-plane vectors, that in turn are used to rotate the velocity vector
    patchNormal_.setSize(nPatches);
    patchTriFaces_.setSize(nPatches);
    sputterdParcels_.setSize(nPatches,Field<scalar>(nSpecies,0.0));
    totalSputterdParcels_.setSize(nPatches,Field<scalar>(nSpecies,0.0));
    runningAverage_.setSize(nPatches,Field<scalar>(nSpecies,0.0));

    forAll(associatedPatches, i)
    {
        label patchId = associatedPatches[i];
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patchId];
        patchTriFaces_[patchId].setSize(patch.size());
        const pointField& points = patch.points();

        DynamicList<face> triFace(2);
        DynamicList<face> tris(5);

        const scalarField magSf(mag(patch.faceAreas()));
        //patchAreas_[patchId] = sum(magSf);

        patchNormal_[patchId] = patch.faceAreas()/magSf;

        forAll(patch, facei)
        {
            const face& f = patch[facei];

            triFace.clear();
            tris.clear();
            f.triangles(points, tris);

            forAll(tris, i)
            {
                triFace.append(tris[i]);
            }

            patchTriFaces_[patchId][facei].transfer(triFace);
        }

    }
}

template<class CloudType>
Foam::SputterEvent<CloudType>::SputterEvent(const SputterEvent<CloudType>& im)
:
    BoundaryEvent<CloudType>(im.owner_),
    sputterSpecies_(im.sputterSpecies_),
    initalEnergies_(im.initalEnergies_),
    sputterProbability_(im.sputterProbability_),
    patchNormal_(im.patchNormal_),
    patchTriFaces_(im.patchTriFaces_),
    printInfo_(im.printInfo_),
    sputterdParcels_(im.sputterdParcels_),
    totalSputterdParcels_(im.totalSputterdParcels_),
    runningAverage_(im.runningAverage_),
    averageStart_(im.averageStart_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SputterEvent<CloudType>::~SputterEvent()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::SputterEvent<CloudType>::info()
{
    CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    if(!printInfo_ || cloud.isInitializing())
        return;

    scalar dt = mesh.time().deltaTValue();
    scalar beta = dt/(mesh.time().value()-averageStart_);

    //Calculate running averages and print info and the rate if particles that are injected
    forAll(this->associatedPatches(),i)
    {
        label patchId = this->associatedPatches()[i];
        const polyPatch& patch = mesh.boundaryMesh()[patchId];


        Pstream::listCombineGather(sputterdParcels_[patchId], plusEqOp<scalar>());
        Pstream::listCombineScatter(sputterdParcels_[patchId]);

        forAll(cloud.typeIdList(),j)
        {
            const typename CloudType::parcelType::constantProperties& cP = cloud.constProps(j);

            totalSputterdParcels_[patchId][j] += sputterdParcels_[patchId][j];

            runningAverage_[patchId][j] = (1.0-beta)*runningAverage_[patchId][j] + beta*sputterdParcels_[patchId][j];

            scalar rate = runningAverage_[patchId][j]/dt;
            Info << "[patch: " << patch.name() << "]" << " Sputter rate (" << cloud.typeIdList()[j] << "): " << rate << " particle/s; I: " << rate * cP.charge() << " A" << endl;

            scalar total = totalSputterdParcels_[patchId][j];
            Info << "[patch: " << patch.name() << "]" << " Total sputter (" << cloud.typeIdList()[j] << "): " << total << " particles" << endl;

            sputterdParcels_[patchId][j] = 0.0;
        }
    }

}

template<class CloudType>
void Foam::SputterEvent<CloudType>::collisionEvent(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{
    CloudType& cloud(this->owner());
    Random& rndGen(cloud.rndGen());
    const fvMesh& mesh(cloud.mesh());

    label patchId = p.patch();
    label typeId = p.typeId();

    if(sputterSpecies_[patchId][typeId].size() > 0)//Do we inject a species after collision with a parcel with "typeId"
    {
        forAll(sputterProbability_[patchId][typeId],i)//Go through all propbilities 
        {
                scalar prop = sputterProbability_[patchId][typeId][i];
                if(rndGen.scalar01() < prop)//Inject?
                {
                    label celli = p.cell();
                    const polyPatch& patch = mesh.boundaryMesh()[patchId];
                    label facei = patch.whichFace(p.face());
                    const point pf(p.position());

                    const vector& fC = patch.faceCentres()[facei];//facecenter
                    const face& tf = patchTriFaces_[patchId][facei][0];//first tri of face
                    vector n = -patchNormal_[patchId][facei];

                    // Position perturbed away from face (into domain)
                    const scalar a = rndGen.scalarAB(0.1, 0.5);
                    const vector& pc = mesh.cellCentres()[celli];
                    const vector d =
                        mag((pf - pc) & patchNormal_[patchId][facei])*patchNormal_[patchId][facei];

                    vector position = pf - a*d;


                    // Wall tangential unit vector. Use the direction between the
                    // face centre and the first vertex in the list
                    vector t1 = fC - (mesh.points()[tf[0]]);
                    t1 /= mag(t1);

                    // Other tangential unit vector.  Rescaling in case face is not
                    // flat and n and t1 aren't perfectly orthogonal
                    vector t2 = n^t1;
                    t2 /= mag(t2);


                    label injectedType = sputterSpecies_[patchId][typeId][i];//TypeId of the parcel to be injected

                    scalar temperature;
                    if(initalEnergies_[patchId][typeId][i] < 0.0)//If energy is negative use the energy of the parcel that hit the patch
                    {
                        //this one includes the dirft velocity, thats an error FIXME...
                        temperature = (p.U()&p.U())*cloud.constProps(typeId).mass()/(3*constant::physicoChemical::k.value());
                    }
                    else
                        temperature = initalEnergies_[patchId][typeId][i]/*eV*/*constant::electromagnetic::e.value()/constant::physicoChemical::k.value();

                    scalar mass = cloud.constProps(injectedType).mass();

                    //3D velocity away from the patch into the domain (uses Maxwell-Boltzmann distribution)
                    scalar R = cloud.rndGen().scalar01();
                    vector U = sqrt(constant::physicoChemical::k.value()*temperature/mass)
                            *(
                                rndGen.scalarNormal()*t1
                                + rndGen.scalarNormal()*t2
                             )
                            + sqrt(constant::physicoChemical::k.value()*temperature/mass)*sqrt(-2.0*log((1.0-R)))*n;

                    //Add the particle
                    typename CloudType::parcelType* newP = cloud.addNewParcel
                    (
                         position,
                         celli,
                         U,
                         injectedType
                    );
                    //Set back the velocity by half a time step using the current time fraction of the particle that hit the patch and setup the weight
                    newP->syncVelocityAtBoundary(cloud, p.stepFraction()-0.5);
                    newP->nParticle() = p.nParticle();

                    sputterdParcels_[patchId][injectedType]+=p.nParticle();

                    cloud.boundaryModels().onEjection(*newP,patchId);
                    cloud.boundaryEventModels().onEjection(*newP,patchId);

                    //Move the new particle only the fraction of this timestep
                    newP->stepFraction() = p.stepFraction();
            }
        }
    }
}


// ************************************************************************* //
