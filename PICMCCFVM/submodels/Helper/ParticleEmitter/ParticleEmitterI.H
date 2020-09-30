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
inline const CloudType& Foam::ParticleEmitter<CloudType>::cloud() const
{
    return cloud_;
}


template<class CloudType>
inline CloudType& Foam::ParticleEmitter<CloudType>::cloud()
{
    return cloud_;
}

template<class CloudType>
inline const Foam::vectorList& Foam::ParticleEmitter<CloudType>::patchNormals() const
{
    return patchNormals_;
}

template<class CloudType>
inline const Foam::label& Foam::ParticleEmitter<CloudType>::patchId() const
{
    return patchId_;
}

template<class CloudType>
inline const Foam::scalar& Foam::ParticleEmitter<CloudType>::patchArea() const
{
    return patchArea_;
}
