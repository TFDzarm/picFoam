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

#include "fvCFD.H"
#include "BorisPusher.H"
#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BorisPusher<CloudType>::BorisPusher
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ParticlePusher<CloudType>(cloud),
    calcMagneticRotation_(false)
{
    //Check the magentic field if we have to calculate rotations
    if(cloud.maxwellSolver().type() == "ElectroStatic" || cloud.maxwellSolver().type() == "none")
    {
        const volVectorField& B(cloud.magneticField());
        forAll(B,i)
        {
            if(mag(B[i]) != 0.0) {
                calcMagneticRotation_ = true;
                break;
            }
        }
       if(!calcMagneticRotation_)
           Info << "|->    No magnetic field, deactivating calculation of the velocity vector rotation" << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BorisPusher<CloudType>::~BorisPusher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BorisPusher<CloudType>::updateVelocity(typename CloudType::parcelType& p, const scalar trackTime)
{
    CloudType& cloud(this->owner());

    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());

    vector& U(p.U());


    const label celli = p.cell();
    const scalar dt2 = 0.5*trackTime;

    const scalar mass = cloud.constProps(p.typeId()).mass();
    const scalar charge = p.charge();

    //update rU
    vector rU = U / Foam::sqrt(1.0-Foam::sqr(Foam::mag(U)/constant::universal::c.value()));
    const vector oldrU = rU;

    const scalar charge_dt2_over_mass = charge*dt2/mass;
    const vector rUsm = charge_dt2_over_mass*E;

    //half step
    rU += rUsm;

    //rotation if |B| is zero everywhere there is no need for these operations
    if(calcMagneticRotation_)
    {
        scalar lF = ( Foam::sqrt(1.0+Foam::sqr((Foam::mag(oldrU)/constant::universal::c.value()))) + Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value()))) ) / 2.0;
        //scalar lF =  Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value()))); //Smilie uses this...
        const vector t = charge_dt2_over_mass*B[celli]/lF;//~150ns
        //vector t = (B[celli]/mag(B[celli]))*Foam::tan(charge*dt*mag(B[celli])/(2*lF*mass));//~500ns div in 10^-14

        const vector Udash = rU + (rU ^ t);
        rU +=  (Udash ^ (2.0*t/(1.0+(t&t))));
    }

    //half step
    rU += rUsm;

    //apply velocity change
    U = rU / Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value())));
}

template<class CloudType>
vector Foam::BorisPusher<CloudType>::correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const
{
    const CloudType& cloud(this->owner());

    //const volVectorField& E(cloud.electricField());
    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());


    const label celli = p.cell();
    const scalar dt2 = 0.5*trackTime;
    const scalar mass = cloud.constProps(p.typeId()).mass();
    const scalar charge = p.charge();

    //update rU
    vector rU = p.U() / Foam::sqrt(1.0-Foam::sqr(Foam::mag(p.U())/constant::universal::c.value()));
    const vector oldrU = rU;

    const scalar charge_dt2_over_mass = charge*dt2/mass;
    const vector rUsm = charge_dt2_over_mass*E;//E[celli];
    //half step
    rU += rUsm;

    //rotation if |B| is zero everywhere there is no need for these operations
    if(calcMagneticRotation_)
    {
        scalar lF = ( Foam::sqrt(1.0+Foam::sqr((Foam::mag(oldrU)/constant::universal::c.value()))) + Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value()))) ) / 2.0;
        const vector t = charge_dt2_over_mass*B[celli]/lF;//~150ns

        const vector Udash = rU + (rU ^ t);
        rU +=  (Udash ^ (2.0*t/(1.0+(t&t))));
    }

    //half step
    rU += rUsm;

    //return velocity change
    return rU / Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value())));
}

// ************************************************************************* //
