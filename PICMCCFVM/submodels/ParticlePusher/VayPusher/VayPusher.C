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
#include "VayPusher.H"
#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VayPusher<CloudType>::VayPusher
(
    const dictionary& dict,
    CloudType& cloud
)
:
    ParticlePusher<CloudType>(cloud)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VayPusher<CloudType>::~VayPusher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VayPusher<CloudType>::updateVelocity(typename CloudType::parcelType& p, const scalar trackTime)
{
    CloudType& cloud(this->owner());

    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());

    vector& U(p.U());

    label celli = p.cell();

    scalar dt2 = 0.5*trackTime;

    scalar mass = cloud.constProps(p.typeId()).mass();
    scalar charge = p.charge();

    scalar lF = 1.0/Foam::sqrt(1.0-Foam::sqr(Foam::mag(p.U())/constant::universal::c.value()));
	//update rU
    vector rU = U*lF;

    vector Ud = rU + charge*trackTime/mass*E;//E[celli]

	vector T = charge*dt2/mass*B[celli];

    Ud += (U^T);

	scalar alpha = 1.0+Foam::sqr(Foam::mag(Ud)/constant::universal::c.value());
	scalar T2 = (T&T);
	
	scalar s = alpha - T2;
	scalar us2 = Foam::pow((Ud&T)/constant::universal::c.value(),2.0);

	alpha = 1.0/Foam::sqrt(0.5*(s + Foam::sqrt(s*s + 4.0*( T2 + us2 ))));
	
	T *= alpha;

	s = 1.0/(1.0+(T&T));
	alpha = Ud&T;

	rU = s*(Ud+(alpha*T)+(Ud^T));//Usm
	
	U = rU / Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value())));

}


template<class CloudType>
vector Foam::VayPusher<CloudType>::correctedVelocity(const typename CloudType::parcelType& p, const scalar trackTime) const
{
    const CloudType& cloud(this->owner());

    const vector E = cloud.eFieldWeighting().getFieldVector(p.coordinates(),p.currentTetIndices());
    const volVectorField& B(cloud.magneticField());

    vector U = p.U();

    label celli = p.cell();

    scalar dt2 = 0.5*trackTime;

    scalar mass = cloud.constProps(p.typeId()).mass();
    scalar charge = p.charge();

    scalar lF = 1.0/Foam::sqrt(1.0-Foam::sqr(Foam::mag(p.U())/constant::universal::c.value()));
    //update rU
    vector rU = U*lF;

    vector Ud = rU + charge*trackTime/mass*E;

    vector T = charge*dt2/mass*B[celli];

    Ud += (U^T);

    scalar alpha = 1.0+Foam::sqr(Foam::mag(Ud)/constant::universal::c.value());
    scalar T2 = (T&T);

    scalar s = alpha - T2;
    scalar us2 = Foam::pow((Ud&T)/constant::universal::c.value(),2.0);

    alpha = 1.0/Foam::sqrt(0.5*(s + Foam::sqrt(s*s + 4.0*( T2 + us2 ))));

    T *= alpha;

    s = 1.0/(1.0+(T&T));
    alpha = Ud&T;

    rU = s*(Ud+(alpha*T)+(Ud^T));//Usm

    return rU / Foam::sqrt(1.0+Foam::sqr((Foam::mag(rU)/constant::universal::c.value())));
}


// ************************************************************************* //
