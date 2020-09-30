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

#include "NoIonization.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoIonization<CloudType>::NoIonization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    IonizationModel<CloudType>(cloud)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoIonization<CloudType>::~NoIonization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NoIonization<CloudType>::initIonization(label pTypeId, label qTypeId, scalar nP, scalar nQ)
{}

template<class CloudType>
void Foam::NoIonization<CloudType>::prepareNumberDensities(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ)
{}

template<class CloudType>
void Foam::NoIonization<CloudType>::prepareIonization(typename CloudType::parcelType& pP, typename CloudType::parcelType& pQ)
{}

template<class CloudType>
bool Foam::NoIonization<CloudType>::ionize()
{
    return false;
}

// ************************************************************************* //
