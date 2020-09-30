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
#include "VariableHardSphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VariableHardSphere<CloudType>::VariableHardSphere
(
    const dictionary& dict,
    CloudType& cloud
)
:
    TotalCrossSectionModel<CloudType>(dict,cloud,typeName),
    Tref_(readScalar(this->coeffDict().lookup("Tref")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VariableHardSphere<CloudType>::~VariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
scalar Foam::VariableHardSphere<CloudType>::sigmaTcR
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

    scalar omegaPQ =
        0.5
       *(
            cloud.constProps(typeIdP).omega()
          + cloud.constProps(typeIdQ).omega()
        );

    scalar cR = mag(pP.U() - pQ.U());

    if (cR < vSmall)
    {
        return 0;
    }

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}

template<class CloudType>
scalar Foam::VariableHardSphere<CloudType>::sigmaTcR
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

    scalar omegaPQ =
        0.5
       *(
            cloud.constProps(typeIdP).omega()
          + cloud.constProps(idQ).omega()
        );

    scalar cR = mag(pP.U() - Uq);

    if (cR < vSmall)
    {
        return 0;
    }

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(idQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}
// ************************************************************************* //

