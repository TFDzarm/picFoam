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

#include "VolumeWeighting.H"
#include "pointConstraints.H"
#include "coupledPointPatchField.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeWeighting<CloudType>::VolumeWeighting
(
    const word& fieldName,
    const dictionary& dict,
    CloudType& cloud
)
:
    ChargeDistribution<CloudType>(fieldName, cloud),
    vertexData_(cloud.mesh().nPoints(),0.0),
    ccData_(cloud.mesh().nCells(),0.0),
    weights_()
{
    const polyMesh& mesh(cloud.mesh());

    Field<scalar> sumWeights(mesh.nPoints(),0.0);

    const pointField& points = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();

    weights_.clear();
    weights_.setSize(mesh.nCells());

    //For every cell calculate the distance from its center to each point constructing the cell and save it as weight
    const labelListList& cellPoints = mesh.cellPoints();
    forAll(cellPoints,celli)
    {
        const labelList& pointList = cellPoints[celli];
        weights_[celli].setSize(pointList.size(),0.0);

        forAll(pointList,idx)
        {
            label pointi = pointList[idx];

            //Distance cell center and vertex
            weights_[celli][idx] = mag(points[pointi] - cellCentres[celli]);
            //Sum of all distances
            sumWeights[pointi] += weights_[celli][idx];
        }
    }

    //Sync the sum, accounting for periodic boundaries
    mesh.globalData().syncPointData
    (
        sumWeights,
        plusEqOp<scalar>(),
        mapDistribute::transform()
    );

    //Calculate the weight
    forAll(cellPoints,celli)
    {
        const labelList& pointList = cellPoints[celli];
        forAll(pointList,idx)
        {
            label pointi = pointList[idx];
            weights_[celli][idx] /= sumWeights[pointi];
        }
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeWeighting<CloudType>::~VolumeWeighting()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::VolumeWeighting<CloudType>::add
(
    const typename CloudType::parcelType& parcel
)
{
    const barycentric& coordinates = parcel.coordinates();
    const tetIndices& tetIs = parcel.currentTetIndices();
    const triFace& triIs(tetIs.faceTriIs(this->cloud().mesh()));

    const scalar value = parcel.charge()*parcel.nParticle();

    //This is the fraction of the total charge in celli that belongs only to celli
    ccData_[tetIs.cell()] += coordinates[0]*value;

    //Add the fraction of the charge to the vertices
    for(label i = 0; i < 3; i ++)
    {
        vertexData_[triIs[i]] += coordinates[i+1]*value;
    }
}

template<class CloudType>
void Foam::VolumeWeighting<CloudType>::reset()
{
    ChargeDistribution<CloudType>::reset();//reset field_

    forAll(vertexData_,i)
       vertexData_[i] = 0.0;
    forAll(ccData_,i)
       ccData_[i] = 0.0;
}

template<class CloudType>
void Foam::VolumeWeighting<CloudType>::update()
{
    //COM the weighted charges, accounts for periodic boundaries
    this->cloud().mesh().globalData().syncPointData
    (
        vertexData_,
        plusEqOp<scalar>(),
        mapDistribute::transform()
    );

    const scalarField& cellVolume = this->cloud().mesh().cellVolumes();
    const labelListList& cellPoints = this->cloud().mesh().cellPoints();

    //Weight back the charge from the vertices to the cell center
    forAll(cellPoints,celli)
    {
        const scalar& V = cellVolume[celli];
        const labelList& pointList = cellPoints[celli];

        this->field_[celli] += ccData_[celli]/V;
        forAll(pointList,idx)
        {
            label pointi = pointList[idx];
            this->field_[celli] += weights_[celli][idx]*vertexData_[pointi]/V;
        }
    }
}

// ************************************************************************* //
