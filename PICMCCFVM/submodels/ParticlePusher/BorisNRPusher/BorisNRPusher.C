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
#include "BorisNRPusher.H"
#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BorisNRPusher<CloudType>::BorisNRPusher
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ParticlePusher<CloudType>(cloud),
    calcMagneticRotation_(false)
{
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
Foam::BorisNRPusher<CloudType>::~BorisNRPusher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BorisNRPusher<CloudType>::updateVelocity(typename CloudType::parcelType& p, const scalar trackTime)
{
    CloudType& cloud(this->owner());

    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());

    vector& U(p.U());

    label celli = p.cell();
    scalar dt2 = 0.5*trackTime;

    scalar mass = cloud.constProps(p.typeId()).mass();
    scalar charge = p.charge();
    scalar alpha = charge*dt2/mass;


    //half step E
    U += E*alpha;

    //rotation: if |B| is zero everywhere there is no need for these operations
    if(calcMagneticRotation_)
    {
        vector t = B[celli]*alpha;//~150ns
        vector Udash = U + (U ^ t);
        U +=  (Udash ^ (2.0*t/(1.0+(t&t))));
    }

    //half step E
    U += E*alpha;
}



template<class CloudType>
vector Foam::BorisNRPusher<CloudType>::correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const
{
    const CloudType& cloud(this->owner());

    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());

    vector U = p.U();

    label celli = p.cell();
    scalar dt2 = 0.5*trackTime;

    scalar mass = cloud.constProps(p.typeId()).mass();
    scalar charge = p.charge();
    scalar alpha = charge*dt2/mass;


    //half step E
    U += E*alpha;

    //rotation: if |B| is zero everywhere there is no need for these operations
    if(calcMagneticRotation_)
    {
        vector t = B[celli]*alpha;//~150ns
        vector Udash = U + (U ^ t);
        U +=  (Udash ^ (2.0*t/(1.0+(t&t))));
    }

    //half step E
    U += E*alpha;
    return U;
}


// ************************************************************************* //
