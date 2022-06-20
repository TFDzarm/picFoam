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

#include "EquipartitionInitialization.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EquipartitionInitialization<CloudType>::EquipartitionInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName),
    cellZoneIndex_(-1)
{
    const fvMesh& mesh(cloud.mesh());
    const dictionary& coeffDict = this->coeffDict();

    word cellZone = coeffDict.lookupOrDefault<word>("cellZone","none");
    if(cellZone != "none") {
        cellZoneIndex_ = mesh.cellZones().findIndex(cellZone);

        if(cellZoneIndex_ < 0)
            FatalErrorInFunction << "Cannot find cellzone \"" << cellZone << "\"" << nl << abort(FatalError);

        Info << "Initializing in cellZone " << cellZone << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EquipartitionInitialization<CloudType>::~EquipartitionInitialization()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::EquipartitionInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());

    const vector velocity(this->coeffDict().lookup("velocity"));

    if(cellZoneIndex() >= 0)
    {
        const labelList& cellLables = cloud.mesh().cellZones()[cellZoneIndex()];
        forAll(cellLables,i)
        {
            label celli = cellLables[i];
            addParticles(celli, velocity);
        }
    }
    else
    {
        forAll(cloud.mesh().cells(), celli)
        {
            addParticles(celli, velocity);
        }
    }
}

template<class CloudType>
void Foam::EquipartitionInitialization<CloudType>::addParticles(label celli, const vector& velocity)
{
    CloudType& cloud(this->cloud());

    List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
    (
        cloud.mesh(),
        celli
    );

    forAll(cellTets, tetI)
    {
        const tetIndices& cellTetIs = cellTets[tetI];
        tetPointRef tet = cellTetIs.tet(cloud.mesh());
        scalar tetVolume = tet.mag();

        forAll(cloud.typeIdList(), typeId)
        {
            const typename CloudType::parcelType::constantProperties& cP =
            cloud.constProps(typeId);

            scalar numberDensity = this->numberDensities_[typeId];

            // Calculate the number of particles required
            scalar particlesRequired = numberDensity*tetVolume / cP.nParticle();

            // Only integer numbers of particles can be inserted
            label nParticlesToInsert = label(particlesRequired);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of particlesRequired
            if
            (
                (particlesRequired - nParticlesToInsert)
              > cloud.rndGen().scalar01()
            )
            {
                nParticlesToInsert++;
            }

            for (label pI = 0; pI < nParticlesToInsert; pI++)
            {
                point p = tet.randomPoint(cloud.rndGen());
                meshTools::constrainDirection(cloud.mesh(), cloud.mesh().solutionD(), p);

                vector U = cloud.equipartitionLinearVelocity
                (
                    this->temperatures_[typeId],
                    cP.mass()
                );

                U += velocity;

                cloud.addNewParcel(p, celli, U, typeId);
            }
        }
    }
}

template<class CloudType>
const Foam::label& Foam::EquipartitionInitialization<CloudType>::cellZoneIndex() const
{
    return cellZoneIndex_;
}

// ************************************************************************* //

