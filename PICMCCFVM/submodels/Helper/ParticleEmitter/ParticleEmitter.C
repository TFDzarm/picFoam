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

#include "ParticleEmitter.H"
#include "Pstream.H"
#include "fvMesh.H"
#include "triPointRef.H"
#include "constants.H"
#include "tetIndices.H"
#include "meshTools.H"
#include "polyMeshTetDecomposition.H"

template<class CloudType>
Foam::ParticleEmitter<CloudType>::ParticleEmitter(CloudType& cloud) :
    cloud_(cloud),
    patchId_(-1),
    patchArea_(0.0),
    patchNormals_(),
    cellOwners_(),
    triFaces_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0),
    getVelocity_(nullptr),
    velocityModel_(VelocityModel::vmUnknown)
{

}

template<class CloudType>
void Foam::ParticleEmitter<CloudType>::initilizeParticleEmitter(label patchId, VelocityModel model)
{
    patchId_ = patchId;
    velocityModel_ = model;

    //Select the chosen model
    switch(model)
    {
        case vmBirdMaxwellianFlux:
        getVelocity_ = &ParticleEmitter::getVelocity0;
        break;
        case vmMostProbableSpeed:
        getVelocity_ = &ParticleEmitter::getVelocity1;
        break;
        case vmHalfMaxwellianFlux:
        getVelocity_ = &ParticleEmitter::getVelocity2;
        break;
        case vmMaxwellianFlux:
        getVelocity_ = &ParticleEmitter::getVelocity3;
        break;
        default:
        FatalErrorInFunction << "Unknown velocity model" << endl;
    }

    const polyPatch& patch = cloud_.mesh().boundaryMesh()[patchId_];
    const pointField& points = patch.points();

    cellOwners_ = patch.faceCells();

    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFaces(2*patch.size());
    DynamicList<face> tris(5);

    triMagSf.append(0.0);

    forAll(patch, facei)
    {
        const face& f = patch[facei];

        tris.clear();
        f.triangles(points, tris);

        forAll(tris, i)
        {
            triToFace.append(facei);
            triFaces.append(tris[i]);
            triMagSf.append(tris[i].mag(points));
        }
    }
    forAll(sumTriMagSf_, i)
    {
        sumTriMagSf_[i] = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); i++)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // Transfer to persistent storage
    triFaces_.transfer(triFaces);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // Convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); i++)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }

    const scalarField magSf(mag(patch.faceAreas()));
    patchArea_ = sum(magSf);
    patchNormals_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());

    Info << "|->    Initilized on patch " << patch.name() << ", area = " << patchArea_ << endl;
}

template<class CloudType>
typename CloudType::parcelType* Foam::ParticleEmitter<CloudType>::emitParticleExplicit(scalar temperature, label typeId, scalar stepFraction, scalar nParticle)
{
    const fvMesh& mesh(cloud_.mesh());
    Random& rndGen(cloud_.rndGen());

    vector position;
    label cellOwner;
    label tetFacei;
    label tetPti;

    //Pick random tri

    label trii = rndGen.sampleAB<label>(0,triToFace_.size());

    // Set cellOwner
    label facei = triToFace_[trii];
    cellOwner = cellOwners_[facei];//use

    // Find random point in triangle
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
    const pointField& points = patch.points();
    const face& tf = triFaces_[trii];
    const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
    const point pf(tri.randomPoint(rndGen));

    // Position perturbed away from face (into domain)
    const scalar a = rndGen.scalarAB(0.1, 0.5);
    const vector& pc = mesh.cellCentres()[cellOwner];
    const vector d =
        mag((pf - pc) & patchNormals_[facei])*patchNormals_[facei];

    position = pf - a*d;//use

    // Try to find tetFacei and tetPti in the current position
    mesh.findTetFacePt(cellOwner, position, tetFacei, tetPti);

    // tetFacei and tetPti not found, check if the cell has changed
    if (tetFacei == -1 ||tetPti == -1)
    {
        mesh.findCellFacePt(position, cellOwner, tetFacei, tetPti);
    }

    // Both searches failed, choose a random position within
    // the original cell
    if (tetFacei == -1 ||tetPti == -1)
    {
        // Reset cellOwner
        cellOwner = cellOwners_[facei];
        const scalarField& V = mesh.V();

        // Construct cell tet indices
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cellOwner);

        // Construct cell tet volume fractions
        scalarList cTetVFrac(cellTetIs.size(), 0.0);
        for (label teti=1; teti<cellTetIs.size()-1; teti++)
        {
            cTetVFrac[teti] =
                cTetVFrac[teti-1]
              + cellTetIs[teti].tet(mesh).mag()/V[cellOwner];
        }
        cTetVFrac.last() = 1;

        // Set new particle position
        const scalar volFrac = rndGen.scalar01();
        label teti = 0;
        forAll(cTetVFrac, vfI)
        {
            if (cTetVFrac[vfI] > volFrac)
            {
                teti = vfI;
                break;
            }
        }

        position = cellTetIs[teti].tet(mesh).randomPoint(rndGen);
        tetFacei = cellTetIs[teti].face();
        tetPti = cellTetIs[teti].tetPt();
    }

    const vector& fC = patch.faceCentres()[facei];

    // Normal unit vector *negative* so normal is pointing into the
    // domain
    vector n = -patchNormals_[facei];

    // Wall tangential unit vector. Use the direction between the
    // face centre and the first vertex in the list
    vector t1 = fC - (mesh.points()[tf[0]]);
    t1 /= mag(t1);

    // Other tangential unit vector.  Rescaling in case face is not
    // flat and n and t1 aren't perfectly orthogonal
    vector t2 = n^t1;
    t2 /= mag(t2);

    //Sample the velocity from the model
    vector U = (this->*getVelocity_)(n, t1, t2, temperature, typeId);

    //Create the particle
    meshTools::constrainDirection(mesh, mesh.solutionD(), position);//2D and 1D constrains
    typename CloudType::parcelType* p = cloud_.addNewParcel
                (
                     position,
                     cellOwner,
                     U,
                     typeId
                );
    if(nParticle > 0.0)
        p->nParticle() = nParticle;

    p->syncVelocityAtBoundary(cloud_, stepFraction-0.5);
    //move new particle only rest of the timestep
    p->stepFraction() = stepFraction;

    //Inform all model of the new particle
    cloud_.boundaryModels().onEjection(*p,patchId_);
    cloud_.boundaryEventModels().onEjection(*p,patchId_);
    return p;

}

template<class CloudType>
typename CloudType::parcelType* Foam::ParticleEmitter<CloudType>::emitParticle(scalar temperature, label typeId, scalar nParticle)
{
    const fvMesh& mesh(cloud_.mesh());

    scalar areaFraction = cloud_.rndGen().globalScalar01()*patchArea_;
    vector position;
    label cellOwner;
    label tetFacei;
    label tetPti;

    if (cellOwners_.size() > 0)
    {
        // Determine which processor to inject from
        label proci = 0;
        forAllReverse(sumTriMagSf_, i)
        {
            if (areaFraction >= sumTriMagSf_[i])
            {
                proci = i;
                break;
            }
        }

        if (Pstream::myProcNo() == proci)
        {
            // Find corresponding decomposed face triangle
            label trii = 0;
            scalar offset = sumTriMagSf_[proci];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    trii = i;
                    break;
                }
            }

            // Set cellOwner
            label facei = triToFace_[trii];
            cellOwner = cellOwners_[facei];//use

            // Find random point in triangle
            const polyPatch& patch = mesh.boundaryMesh()[patchId_];
            const pointField& points = patch.points();
            const face& tf = triFaces_[trii];
            const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
            const point pf(tri.randomPoint(cloud_.rndGen()));

            // Position perturbed away from face (into domain)
            const scalar a = cloud_.rndGen().scalarAB(0.1, 0.5);
            const vector& pc = mesh.cellCentres()[cellOwner];
            const vector d =
                mag((pf - pc) & patchNormals_[facei])*patchNormals_[facei];

            position = pf - a*d;//use

            // Try to find tetFacei and tetPti in the current position
            mesh.findTetFacePt(cellOwner, position, tetFacei, tetPti);

            // tetFacei and tetPti not found, check if the cell has changed
            if (tetFacei == -1 ||tetPti == -1)
            {
                mesh.findCellFacePt(position, cellOwner, tetFacei, tetPti);
            }

            // Both searches failed, choose a random position within
            // the original cell
            if (tetFacei == -1 ||tetPti == -1)
            {
                // Reset cellOwner
                cellOwner = cellOwners_[facei];
                const scalarField& V = mesh.V();

                // Construct cell tet indices
                const List<tetIndices> cellTetIs =
                    polyMeshTetDecomposition::cellTetIndices(mesh, cellOwner);

                // Construct cell tet volume fractions
                scalarList cTetVFrac(cellTetIs.size(), 0.0);
                for (label teti=1; teti<cellTetIs.size()-1; teti++)
                {
                    cTetVFrac[teti] =
                        cTetVFrac[teti-1]
                      + cellTetIs[teti].tet(mesh).mag()/V[cellOwner];
                }
                cTetVFrac.last() = 1;

                // Set new particle position
                const scalar volFrac = cloud_.rndGen().scalar01();
                label teti = 0;
                forAll(cTetVFrac, vfI)
                {
                    if (cTetVFrac[vfI] > volFrac)
                    {
                        teti = vfI;
                        break;
                    }
                }
                position = cellTetIs[teti].tet(mesh).randomPoint(cloud_.rndGen());
                tetFacei = cellTetIs[teti].face();
                tetPti = cellTetIs[teti].tetPt();
            }

            const vector& fC = patch.faceCentres()[facei];

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = -patchNormals_[facei];

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[tf[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            //Sample the velocity from the model
            vector U = (this->*getVelocity_)(n, t1, t2, temperature, typeId);


            //Create the particle
            meshTools::constrainDirection(mesh, mesh.solutionD(), position);//2D and 1D constrains
            typename CloudType::parcelType* p = cloud_.addNewParcel
                        (
                             position,
                             cellOwner,
                             U,
                             typeId
                        );
            if(nParticle > 0.0)
                p->nParticle() = nParticle;

            //Sync Particle if requierd and update boundary models!
            p->syncVelocityAtBoundary(cloud_, -0.5);

            //Inform all models about the new particle
            cloud_.boundaryModels().onEjection(*p,patchId_);
            cloud_.boundaryEventModels().onEjection(*p,patchId_);
            return p;
          }
    }
    return nullptr;//not all procs return a particle
}

template<class CloudType>
Foam::vector Foam::ParticleEmitter<CloudType>::getVelocity(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId)
{
    return (this->*getVelocity_)(normal, tangent1, tangent2, temperature, typeId);//Call the function
}

//BirdMaxwellianFlux
template<class CloudType>
Foam::vector Foam::ParticleEmitter<CloudType>::getVelocity0(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId)
{
    scalar mass = cloud_.constProps(typeId).mass();

    vector rndOrthogonalU = sqrt(constant::physicoChemical::k.value()*temperature/mass)
            *(
                 cloud_.rndGen().scalarNormal()*tangent1
               + cloud_.rndGen().scalarNormal()*tangent2
             );

    scalar mostProbableSpeed
    (
        cloud_.maxwellianMostProbableSpeed
        (
            temperature,
            mass
        )
    );

    scalar P = -1;
    // Normalised candidates for the normal direction velocity
    // component
    scalar uNormal;

    // Select a velocity using Bird eqn 12.5
    do
    {
        uNormal = 3.0*(2.0*cloud_.rndGen().scalar01() - 1.0);


        if (uNormal < 0.0)
            P = -1;
        else
            P = 2.0*uNormal/sqrt(2.0)*exp(0.5 - sqr(uNormal));

    } while (P < cloud_.rndGen().scalar01());

    vector U = rndOrthogonalU + mostProbableSpeed*uNormal*normal;

    return U;
}

//MostProbableSpeed
template<class CloudType>
Foam::vector Foam::ParticleEmitter<CloudType>::getVelocity1(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId)
{
    scalar mass = cloud_.constProps(typeId).mass();

    scalar mostProbableSpeed
    (
        cloud_.maxwellianMostProbableSpeed
        (
            temperature,
            mass
        )
    );

    vector U = mostProbableSpeed*normal;

    return U;
}

//HalfMaxwellianFlux
template<class CloudType>
Foam::vector Foam::ParticleEmitter<CloudType>::getVelocity2(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId)
{
    scalar mass = cloud_.constProps(typeId).mass();

    vector rndOrthogonalU = sqrt(constant::physicoChemical::k.value()*temperature/mass)
            *(
                 cloud_.rndGen().scalarNormal()*tangent1
               + cloud_.rndGen().scalarNormal()*tangent2
             );

    //magnitude of normal distribution [-inv;inv] => [0;inv] is twice as likely so we got a half-maxwellian distribution in x
    vector U = rndOrthogonalU + Foam::mag(cloud_.rndGen().scalarNormal())*sqrt(constant::physicoChemical::k.value()*temperature/mass)*normal;

    return U;

}

//MaxwellianFlux
template<class CloudType>
Foam::vector Foam::ParticleEmitter<CloudType>::getVelocity3(const vector& normal, const vector& tangent1, const vector& tangent2, scalar temperature, label typeId)
{
    namespace cu = constant::universal;

    scalar mass = cloud_.constProps(typeId).mass();

    //thermal velocity
    scalar vt = sqrt(constant::physicoChemical::k.value()*temperature/mass);

    vector rndOrthogonalU = vt
            *(
                 cloud_.rndGen().scalarNormal()*tangent1
               + cloud_.rndGen().scalarNormal()*tangent2
             );
    scalar R = cloud_.rndGen().scalar01();

    //cutoffs -> here we use no cutoffs
    /*
    scalar vL = 0.0;
    scalar vU = cu::c.value();
    scalar v0 = 0.0;
    if(temperature > 0.0)
        vU /= vt;
    scalar expvU2 = ::exp(-0.5*sqr(vU));
    scalar expvL2 = ::exp(-0.5*sqr(vL));

    scalar uNorm = mag(v0) + vt*sqrt(-2.0*log((1.0-R)*expvL2 + R*expvU2));
    */

    //inverse of the maxwellian distribution [0;inv]
    scalar uNorm = vt*sqrt(-2.0*log((1.0-R)));
    return rndOrthogonalU + uNorm*normal;
}
