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

template<class CloudType>
inline const bool& Foam::ElectronNeutralCollisionModel<CloudType>::deleteNeutrals() const
{
    return deleteNeutral_;
}

template<class CloudType>
inline const bool& Foam::ElectronNeutralCollisionModel<CloudType>::updateNeutralVelocity() const
{
    return updateNeutralVelocity_;
}

template<class CloudType>
inline const bool& Foam::ElectronNeutralCollisionModel<CloudType>::createIon() const
{
    return createIon_;
}

template<class CloudType>
inline const bool& Foam::ElectronNeutralCollisionModel<CloudType>::createElectron() const
{
    return createElectron_;
}

template<class CloudType>
inline const CrossSectionList<CloudType,Foam::crossSectionType::ElectronElasticCS>&
    Foam::ElectronNeutralCollisionModel<CloudType>::elasticCrossSections() const
{
    return elasticCSModels_;
}

template<class CloudType>
inline const CrossSectionList<CloudType,Foam::crossSectionType::ElectronExciationCS>&
    Foam::ElectronNeutralCollisionModel<CloudType>::excitationCrossSections() const
{
    return excitationCSModels_;
}

template<class CloudType>
inline const CrossSectionList<CloudType,Foam::crossSectionType::ElectronIonizationCS>&
    Foam::ElectronNeutralCollisionModel<CloudType>::ionizationCrossSections() const
{
    return ionizationCSModels_;
}
