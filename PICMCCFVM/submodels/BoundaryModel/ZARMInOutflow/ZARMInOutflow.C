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

#include "ZARMInOutflow.H"
#include "constants.H"
#include "triPointRef.H"
#include "tetIndices.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ZARMInOutflow<CloudType>::ZARMInOutflow
(
    const dictionary& dict,
    CloudType& cloud,
    const List<label>& associatedPatches
)
:
    BoundaryModel<CloudType>(dict, cloud, typeName, associatedPatches),
    patches_(),
    moleculeTypeIds_(),
    numberDensities_(),
    particleFluxAccumulators_(),
    interactionList_(cloud.mesh().boundaryMesh().size(),pTNoInteraction),
    boundaryT_(),
    boundaryU_()
{
    // Identify which patches to use
    DynamicList<label> patches;

    forAll(associatedPatches, i)
    {
        label p = associatedPatches[i];//patchID

        const polyPatch& patch = cloud.mesh().boundaryMesh()[p];

        if (isType<polyPatch>(patch))
        {
            patches.append(p);

            //Lookup type
            const dictionary& patchPropertiesDict = this->coeffDict().subDict(patch.name());
            const word type(patchPropertiesDict.lookup("type"));
            if(type == "inlet")
                interactionList_[p] = pTInlet;
            else if(type == "inletOutlet")
                interactionList_[p] = pTInletOutlet;
            else
            {
                FatalErrorInFunction
                        << "patch type " << type << " for patch " << patch.name() << " not defined. valid options: [inlet/inletOutlet]" << nl
                        << abort(FatalError);
            }
        }
    }

    patches_.transfer(patches);

    List<wordList> molecules(patches_.size());
    forAll(patches_, p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patches_[p]];
        const patchType& pT = interactionList_[patches_[p]];

        //On the inlet the faceFlux [1/(m^2 s)] is given...
        const dictionary& numberDensitiesDict = (pT == pTInlet) ? this->coeffDict().subDict(patch.name()).subDict("faceFlux") : this->coeffDict().subDict(patch.name()).subDict("numberDensities");
        molecules[p] = wordList(
                    numberDensitiesDict.toc()
                    );
    }

    // Initialise the particleFluxAccumulators_
    particleFluxAccumulators_.setSize(patches_.size());

    forAll(patches_, p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patches_[p]];

        particleFluxAccumulators_[p] = List<Field<scalar>>
        (
            molecules[p].size(),
            Field<scalar>(patch.size(), 0.0)
        );
    }

    moleculeTypeIds_.setSize(patches_.size());
    numberDensities_.setSize(patches_.size());
    boundaryT_.setSize(patches_.size());
    boundaryU_.setSize(patches_.size());

    forAll(patches_, p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patches_[p]];
        const patchType& pT = interactionList_[patches_[p]];

        moleculeTypeIds_[p].setSize(molecules[p].size());
        numberDensities_[p].setSize(molecules[p].size());
        //On the inlet the faceFlux [1/(m^2 s)] is given...
        const dictionary& numberDensitiesDict = (pT == pTInlet) ? this->coeffDict().subDict(patch.name()).subDict("faceFlux") : this->coeffDict().subDict(patch.name()).subDict("numberDensities");

        const dictionary& temperatureDict
        (
            this->coeffDict().subDict(patch.name()).subDict("temperature")
        );

        const dictionary& velocityDict
        (
            this->coeffDict().subDict(patch.name()).subDict("velocity")
        );

        forAll(molecules[p], i)
        {
            numberDensities_[p][i] = readScalar
            (
                numberDensitiesDict.lookup(molecules[p][i])
            );

            //This is one fixed temperature, do we require a field?
            boundaryT_[p][i] = readScalar
            (
                temperatureDict.lookup(molecules[p][i])
            );

            if (boundaryT_[p][i] < small)
            {
                FatalErrorInFunction
                    << "Zero boundary temperature detected, check value for species "
                    << molecules[i] << " at patch " << patch.name() << nl
                    << nl << abort(FatalError);
            }

            boundaryU_[p][i] = velocityDict.lookup(molecules[p][i]);

            label typeId = findIndex(cloud.typeIdList(), molecules[p][i]);
            moleculeTypeIds_[p][i] = typeId;

            if (moleculeTypeIds_[p][i] == -1)
            {
                FatalErrorInFunction
                    << "typeId " << molecules[p][i] << " not defined in cloud." << nl
                    << abort(FatalError);
            }
            //Number density as parcels so we do not need to divide by nParticle each time
            numberDensities_[p][i] /= cloud.constProps(typeId).nParticle();
        }
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ZARMInOutflow<CloudType>::~ZARMInOutflow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ZARMInOutflow<CloudType>::injection()
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    const scalar deltaT = mesh.time().deltaTValue();

    Random& rndGen(cloud.rndGen());

    scalar sqrtPi = sqrt(pi);

    label particlesInserted = 0;

    forAll(patches_, p)
    {
        label patchi = patches_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Add mass to the accumulators.  negative face area dotted with the
        // velocity to point flux into the domain.

        // Take a reference to the particleFluxAccumulator for this patch
        List<Field<scalar>>& pFA = particleFluxAccumulators_[p];

        forAll(pFA, i)
        {
            label typeId = moleculeTypeIds_[p][i];

            scalar mass = cloud.constProps(typeId).mass();

            scalarField mostProbableSpeed
            (
                cloud.maxwellianMostProbableSpeed
                (
                    boundaryT_[patchi][i],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalarField sCosTheta
            (
                (boundaryU_[patchi] & -patch.faceAreas()/mag(patch.faceAreas()))
              / mostProbableSpeed
            );

            //Calculate different pFA's for inlet and outlet patches
            if(interactionList_[p] == pTInlet)
            {
                pFA[i] +=
                    mag(patch.faceAreas())*numberDensities_[p][i]*deltaT;//numberDensities_ == faceFlux 1/(m^2 s)
            }
            else
            {
                pFA[i] +=
                    mag(patch.faceAreas())*numberDensities_[p][i]*deltaT
                    *mostProbableSpeed
                    *(
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                     )
                    /(2.0*sqrtPi);
            }
        }

        forAll(patch, pFI)
        {
            // Loop over all faces as the outer loop to avoid calculating
            // geometrical properties multiple times.

            const face& f = patch[pFI];

            label globalFaceIndex = pFI + patch.start();

            label celli = mesh.faceOwner()[globalFaceIndex];

            const vector& fC = patch.faceCentres()[pFI];

            scalar fA = mag(patch.faceAreas()[pFI]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh,
                globalFaceIndex,
                celli
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                  + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = patch.faceAreas()[pFI];
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[f[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            forAll(pFA, i)
            {
                scalar faceTemperature = boundaryT_[patchi][i];
                const vector& faceVelocity = boundaryU_[patchi][i];

                scalar& faceAccumulator = pFA[i][pFI];

                // Number of whole particles to insert
                label nI = max(label(faceAccumulator), 0);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of faceAccumulator
                if ((faceAccumulator - nI) > rndGen.scalar01())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                label typeId = moleculeTypeIds_[p][i];

                scalar mass = cloud.constProps(typeId).mass();

                for (label i = 0; i < nI; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.scalar01();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle.

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point p = faceTetIs.faceTri(mesh).randomPoint(rndGen);

                    // Velocity generation

                    scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;

                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.scalar01() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.scalar01());

                    vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                       *(
                            rndGen.scalarNormal()*t1
                          + rndGen.scalarNormal()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    //Create a new parcel, sync with the leap frog scheme and inform boundaries of this creation
                    typename CloudType::parcelType* parcel = cloud.addNewParcel
                    (
                        p,
                        celli,
                        U,
                        typeId
                    );
                    parcel->syncVelocityAtBoundary(cloud, -0.5);
                    particleEjection(*parcel,patchi);//directly call
                    cloud.boundaryEventModels().onEjection(*parcel,patchi);

                    particlesInserted++;
                }
            }
        }
    }

    reduce(particlesInserted, sumOp<label>());

    Info<< "[" << typeName << "] Particles inserted              = "
        << particlesInserted << endl;

}

template<class CloudType>
bool Foam::ZARMInOutflow<CloudType>::particleBC(typename CloudType::parcelType& p, typename CloudType::parcelType::trackingData& td)
{ 
    if(interactionList_[p.patch()] == pTInlet)
    {
        //Do a simple specular reflection,
        //this patch is open => there is no actual wall
        p.specularReflect();
        td.requireResync() = true;//Is this correct? If we do not do this charged particles are in sync going forward...
        return true;
    }

    //All particles reaching this patch will be deleted, set this here so BoundaryEvent models know about this
    td.keepParticle = false;
    return false;
}


// ************************************************************************* //
