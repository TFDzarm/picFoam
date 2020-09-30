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

template<class CloudType>
inline const bool& Foam::CoulombCollisionModel<CloudType>::allowIntraCollision() const
{
    return allowIntraCollision_;
}

template<class CloudType>
inline const bool& Foam::CoulombCollisionModel<CloudType>::calculateDebyeLength() const
{
    return calculateDebyeLength_;
}

template<class CloudType>
inline const Foam::PairingAlgorithm<CloudType>&
Foam::CoulombCollisionModel<CloudType>::pairingAlgorithm() const
{
    return pairingAlgorithm_;
}

template<class CloudType>
inline const Foam::IonizationModel<CloudType>&
Foam::CoulombCollisionModel<CloudType>::ionizationModel() const {
    return ionizationModel_;
}

template<class CloudType>
inline Foam::IonizationModel<CloudType>&
Foam::CoulombCollisionModel<CloudType>::ionizationModel() {
    return ionizationModel_();
}

template<class CloudType>
inline const Foam::WeightCorrectionModel<CloudType>&
Foam::CoulombCollisionModel<CloudType>::weightCorrection() const {
    return weightCorrection_;
}

template<class CloudType>
inline Foam::WeightCorrectionModel<CloudType>&
Foam::CoulombCollisionModel<CloudType>::weightCorrection() {
    return weightCorrection_();
}
