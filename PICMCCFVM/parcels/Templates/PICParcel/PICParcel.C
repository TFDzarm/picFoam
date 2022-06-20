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

#include "PICParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
Foam::PICParcel<ParcelType>::move called by the cloud base class

Update the velocity and move the parcel
*/
template<class ParcelType>
template<class TrackCloudType>
bool Foam::PICParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    //reset td
    td.switchProcessor = false;
    td.keepParticle = true;

    //should the parcel move?
    const constantProperties& constProps(cloud.constProps(typeId_));
    if(!constProps.solveMovement())
        return true;

    //Update velocity only once! ::move is called again after processor switch
    if(charge() != 0.0 && this->stepFraction() == 0.0)
    {
        cloud.particlePusher().updateVelocity(*this, trackTime);//call into the particle pusher submodel
    }

    //stop simulation if it becomes unphysical
    if(mag(U_) >= constant::universal::c.value())
    {
        FatalError << "Particle velocity greater than speed of light (" << constant::universal::c.value() << ")! Aborting..." << nl
                   << "typeId: " << typeId_ << nl << "Velocity: " << U_ << " mag(" << mag(U_) << ")" << nl
                   << "|v-c|: " << mag(U_) - constant::universal::c.value() << abort(FatalError);
    }

    //the actual move
    moveForward(cloud, td, trackTime);

    return td.keepParticle;//particle will be removed in the base class if false
}

/*
Foam::PICParcel<ParcelType>::moveForward called by Foam::PICParcel<ParcelType>::move

Move the parcel from cell face to cell face and perform boundary interactions
*/
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::moveForward
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyMesh& mesh = cloud.pMesh();

    // Get cell length scale
    const scalarField& cellLengthScale = cloud.cellLengthScale();


    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    // Total displacement over the time-step
    const vector s = trackTime*Utracking;

    // Trajectory warning
    //Notice! Prints only first particle, there may be faster ones... -> Save max and compare against avg cell size?
    const scalar& l = cellLengthScale[p.cell()]*cloud.warnCellTrajectory();
    if(cloud.printVelocityWarning() && mag(s) > l)
    {
        //FIXME: this is just the current error maybe print an average over all and a max...
        Pout << "    WARNING parcel trajectory (" << mag(s)/cellLengthScale[p.cell()] << ") is larger than the given limit (" << cloud.warnCellTrajectory() << ")" << nl << endl;
        cloud.printVelocityWarning() = false;
    }

    label resyncCount = 0;
    while (td.keepParticle && !td.switchProcessor && p.stepFraction() < 1)
    {
        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        const scalar f = 1 - p.stepFraction();

        p.trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);//move to the next face and call the boundary interaction

        if(td.requireResync() && td.keepParticle)
        {
            resyncCount++;
            if(resyncCount > 1) { //Rarely happens... but the resync can get the parcel stuck, since the parcel should move less then a cell length, more than one sync indicates a stuck particle
                td.requireResync() = false;
                continue;
            }
            //Sync velocity back after reflection
            syncVelocityAtBoundary(cloud, 0.5-p.stepFraction());
            td.requireResync() = false;
        }
    }
}

/*
Foam::PICParcel<ParcelType>::hitPatch called by the base class

If we hit a patch first sync the particle then interact with it
*/
template<class ParcelType>
template<class TrackCloudType>
bool Foam::PICParcel<ParcelType>::hitPatch(TrackCloudType& cloud, trackingData& td)
{
    typename TrackCloudType::parcelType& p = static_cast<typename TrackCloudType::parcelType&>(*this);
    const polyPatch& patch = cloud.mesh().boundaryMesh()[this->patch()];

    //Only do sync on simple polyPatches or wallPolyPatches
    //Alternative: Move hitSymmetryPatch, cylic, ... and require resync but that would be extra work not needed
    if(isType<polyPatch>(patch) || isA<wallPolyPatch>(patch))
        syncVelocityAtBoundary(cloud, p.stepFraction()-0.5);

    //set boundaryVelocity
    td.boundaryVelocity() = U_;

    //Particle boundary condition on patches
    bool hitPatch = cloud.boundaryModels().particleBC(*this, td);


    //Moved this from particleTemplates.C void Foam::particle::hitFace(...) to this place so we can do secondary emission after the hit!
    //If we haven't hit a patch yet and the patch is a wall perform a wall interaction
    if (!hitPatch && isA<wallPolyPatch>(patch))
    {
        wallReflection(cloud, td);
        hitPatch = true;
        //Wall reflection: resync required!
        td.requireResync() = true;
    }

    //Notify event models
    cloud.boundaryEventModels().onCollision(*this, td);

    return hitPatch;
}

/*
Foam::PICParcel<ParcelType>::hitProcessorPatch called by the base class

Set td so the base class will perform a parallel communication
*/
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}
/*
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::hitSymmetryPatch(TrackCloudType&, trackingData& td)
{
    td.requireResync() = true;
    const vector nf = this->normal();
    transformProperties(I - 2.0*nf*nf);
    //ParcelType::hitSymmetryPatch(TrackCloudType& cloud, trackingData& td);
}*/


/*
Foam::PICParcel<ParcelType>::syncVelocityAtBoundary

If syncing is enabled, sync the velocity according to the current move fraction
*/
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::syncVelocityAtBoundary(TrackCloudType& cloud, scalar fraction)
{
    if(cloud.syncVelocityAtBoundary() && charge_ != 0.0)//Syncing makes only sense for charged species, thus we save some instructions if we check the charge.
    {
        scalar syncDT = fraction*cloud.mesh().time().deltaTValue();
        cloud.particlePusher().updateVelocity(*this, syncDT);
    }
}

/*
Foam::PICParcel<ParcelType>::wallReflection

Update fields and call the appropriate wall interaction model
*/
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::wallReflection
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    const label fieldIndex = cloud.fieldCalculation()[typeId_];

    const label wppIndex = this->patch();

    const wallPolyPatch& wpp =
        static_cast<const wallPolyPatch&>
        (
            this->mesh().boundaryMesh()[wppIndex]
        );

    const label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = cloud.pMesh().time().deltaTValue();

    const constantProperties& constProps(cloud.constProps(typeId_));

    scalar m = constProps.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = U_ & nw;

    vector Ut = U_ - U_dot_nw*nw;

    scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, vSmall);

    cloud.rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    cloud.rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    cloud.linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    cloud.momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    if(fieldIndex >= 0)
    {
        cloud.rhoMSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*invMagUnfA;
        cloud.rhoNSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += invMagUnfA;
        cloud.momentumSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
        cloud.linearKESpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += 0.5*m*(U_ & U_)*invMagUnfA;
    }

    // pre-interaction energy
    scalar preIE = 0.5*m*(U_ & U_);

    // pre-interaction momentum
    vector preIMom = m*U_;

    // -------- WallInteraction/Reflection --------
    cloud.wallReflection().correct(*this);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    U_dot_nw = U_ & nw;

    Ut = U_ - U_dot_nw*nw;

    invMagUnfA = 1/max(mag(U_dot_nw)*fA, vSmall);

    cloud.rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    cloud.rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    cloud.linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    cloud.momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    if(fieldIndex >= 0)
    {
        cloud.rhoMSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*invMagUnfA;
        cloud.rhoNSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += invMagUnfA;
        cloud.momentumSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
        cloud.linearKESpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += 0.5*m*(U_ & U_)*invMagUnfA;
    }

    // post-interaction energy
    scalar postIE = 0.5*m*(U_ & U_);

    // post-interaction momentum
    vector postIMom = m*U_;

    scalar deltaQ = nParticle_*(preIE - postIE)/(deltaT*fA);

    vector deltaFD = nParticle_*(preIMom - postIMom)/(deltaT*fA);

    cloud.qBF()[wppIndex][wppLocalFace] += deltaQ;

    cloud.fDBF()[wppIndex][wppLocalFace] += deltaFD;

    //requireResync!
    td.requireResync() = true;
}

/*
Foam::PICParcel<ParcelType>::specularReflect

Simple specular reflection used by some models which are no wall and therefore do not update fields
*/
template<class ParcelType>
void Foam::PICParcel<ParcelType>::specularReflect()
{
    const vector nw = this->normal();

    scalar U_dot_nw = U_ & nw;

    if (U_dot_nw > 0.0)
    {
        U_ -= 2.0*U_dot_nw*nw;
    }
}

/*
Foam::PICParcel<ParcelType>::wallAbsorption

Update fields and remove the particle
*/
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::wallAbsorption
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    vector U = td.boundaryVelocity();
    const label fieldIndex = cloud.fieldCalculation()[typeId_];

    const label wppIndex = this->patch();

    const wallPolyPatch& wpp =
        static_cast<const wallPolyPatch&>
        (
            this->mesh().boundaryMesh()[wppIndex]
        );

    const label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = cloud.pMesh().time().deltaTValue();

    const constantProperties& constProps(cloud.constProps(typeId_));

    scalar m = constProps.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = U & nw;

    vector Ut = U - U_dot_nw*nw;

    scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, vSmall);

    cloud.rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    cloud.rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    cloud.linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U & U)*invMagUnfA;

    cloud.momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    if(fieldIndex >= 0)
    {
        cloud.rhoMSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*invMagUnfA;
        cloud.rhoNSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += invMagUnfA;
        cloud.momentumSpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
        cloud.linearKESpeciesBF(fieldIndex)[wppIndex][wppLocalFace] += 0.5*m*(U & U)*invMagUnfA;
    }

    scalar deltaQ = nParticle_*(0.5*m*(U_ & U_))/(deltaT*fA);
    vector deltaFD = nParticle_*(m*U_)/(deltaT*fA);
    cloud.qBF()[wppIndex][wppLocalFace] += deltaQ;
    cloud.fDBF()[wppIndex][wppLocalFace] += deltaFD;
    //cloud.jBF()[wppIndex][wppLocalFace] += nParticle_*charge_*U_/fA;

    // -------- WallAbsorption --------
    //U = zero;

    td.keepParticle = false;// Deletion
}

//Is not called, overwritten here just to make it clear!
template<class ParcelType>
template<class TrackCloudType>
void Foam::PICParcel<ParcelType>::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{} //We handle wall interaction in ::hitPatch!



template<class ParcelType>
void Foam::PICParcel<ParcelType>::transformProperties(const transformer& transform)
{
    ParcelType::transformProperties(transform);
    U_ = transform.transform(U_);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "PICParcelIO.C"

// ************************************************************************* //
