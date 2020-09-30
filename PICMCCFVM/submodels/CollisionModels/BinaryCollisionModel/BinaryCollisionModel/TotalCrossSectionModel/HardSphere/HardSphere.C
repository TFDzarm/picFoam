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
#include "HardSphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HardSphere<CloudType>::HardSphere
(
    const dictionary& dict,
    CloudType& cloud
)
:
    TotalCrossSectionModel<CloudType>(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HardSphere<CloudType>::~HardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
scalar Foam::HardSphere<CloudType>::sigmaTcR
(
    const typename CloudType::parcelType& pP,
    const typename CloudType::parcelType& pQ
) const
{
    const CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    scalar dPQ =
        0.5
       *(
            cloud.constProps(typeIdP).d()
          + cloud.constProps(typeIdQ).d()
        );

    scalar cR = mag(pP.U() - pQ.U());

    if (cR < vSmall)
    {
        return 0;
    }

    // calculating cross section from Bird, eq. 2.31
    scalar sigmaTPQ =
        pi*dPQ*dPQ;

    return sigmaTPQ*cR;
}

template<class CloudType>
scalar Foam::HardSphere<CloudType>::sigmaTcR
(
    const typename CloudType::parcelType& pP,
    const vector& Uq,
    label idQ
) const
{
    const CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();

    scalar dPQ =
        0.5
       *(
            cloud.constProps(typeIdP).d()
          + cloud.constProps(idQ).d()
        );

    scalar cR = mag(pP.U() - Uq);

    if (cR < vSmall)
    {
        return 0;
    }

    // calculating cross section from Bird, eq. 2.31
    scalar sigmaTPQ =
        pi*dPQ*dPQ;

    return sigmaTPQ*cR;
}
// ************************************************************************* //

