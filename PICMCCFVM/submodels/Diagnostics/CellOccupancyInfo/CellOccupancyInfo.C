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

#include "CellOccupancyInfo.H"
#include "constants.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellOccupancyInfo<CloudType>::CellOccupancyInfo
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DiagnosticInfo<CloudType>(dict,cloud,typeName),
    nCells_(cloud.mesh().nCells())
{
    reduce(nCells_, sumOp<label>());//Total number of cells is fixed, dynamic meshes are not supported at the moment...
}

template<class CloudType>
Foam::CellOccupancyInfo<CloudType>::CellOccupancyInfo(const CellOccupancyInfo<CloudType>& im)
:
    DiagnosticInfo<CloudType>(im.dict_, im.owner_, im.typeName),
    nCells_(im.nCells_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellOccupancyInfo<CloudType>::~CellOccupancyInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CellOccupancyInfo<CloudType>::gatherDiagnostic(const typename CloudType::parcelType& p)
{}

//- Print info
template<class CloudType>
void Foam::CellOccupancyInfo<CloudType>::info()
{
    const CloudType& cloud(this->owner());

    const scalarField& picRhoN = cloud.picRhoN().primitiveField();

    scalar gAverage = 0.0;
    scalar gMin = numeric_limits<scalar>::max();
    scalar gMax = 0.0;//There won't be negativ numbers

    //Loop so we can check multiple stats at once
    forAll(picRhoN, i)
    {
        const scalar& N = picRhoN[i];
        if(N < gMin)
            gMin = N;
        if(N > gMax)
            gMax = N;
        gAverage += N;
    }

    reduce(gAverage, sumOp<scalar>());
    reduce(gMin, minOp<scalar>());
    reduce(gMax, maxOp<scalar>());


    Info << "    Cell occupancy" << nl
         << "    average = " << gAverage/nCells_ << " | maximum = " << gMax << " | minimum = " << gMin << endl;

}

// ************************************************************************* //
