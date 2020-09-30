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

#include "DiagnosticInfo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiagnosticInfo<CloudType>::DiagnosticInfo
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subOrEmptyDict(type + "Coeffs")),
    diagnosticControl_(owner.db().time(),coeffDict_,"diagnostic")
    //writeControl_(owner.db().time(),coeffDict_,"write")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DiagnosticInfo<CloudType>::~DiagnosticInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::DiagnosticInfo<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::DiagnosticInfo<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::DiagnosticInfo<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::DiagnosticInfo<CloudType>::coeffDict() const
{
    return coeffDict_;
}

template<class CloudType>
bool Foam::DiagnosticInfo<CloudType>::shouldExecute()
{
    return diagnosticControl_.execute();
}
/*
template<class CloudType>
bool Foam::DiagnosticInfo<CloudType>::shouldWrite()
{
    return writeControl_.execute();
}*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DiagnosticInfoNew.C"

// ************************************************************************* //
