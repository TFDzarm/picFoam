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

#include "MaxwellianInitialization.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MaxwellianInitialization<CloudType>::MaxwellianInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName),
    driftVelocity_(cloud.typeIdList().size(),Zero),
    thermalVelocity_(cloud.typeIdList().size(),Zero)
{
    forAllConstIter(IDLList<entry>, this->coeffDict(), iter)
    {
        if(iter().isDict())
        {
            label typeId = findIndex(cloud.typeIdList(),iter().keyword());
            if(typeId == -1)
                FatalErrorInFunction << "Undefined typeId " << iter().keyword() << abort(FatalError);

            const dictionary& subDict = iter().dict();
            driftVelocity_[typeId] = subDict.lookup("driftVelocity");
            thermalVelocity_[typeId] = subDict.lookup("thermalVelocity");
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MaxwellianInitialization<CloudType>::~MaxwellianInitialization()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::MaxwellianInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    Random& rndGen(cloud.rndGen());

    forAll(cloud.mesh().cells(), celli)
    {
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
                if((particlesRequired - nParticlesToInsert) > rndGen.scalar01())
                    nParticlesToInsert++;

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(cloud.rndGen());
                    meshTools::constrainDirection(cloud.mesh(), cloud.mesh().solutionD(), p);

                    vector U = thermalVelocity_[typeId];
                    vector n = rndGen.sampleNormal<vector>();
                    U[0] *= n[0];
                    U[1] *= n[1];
                    U[2] *= n[2];

                    U += driftVelocity_[typeId];

                    cloud.addNewParcel(p, celli, U, typeId);
                }

            }
        }
    }
}

// ************************************************************************* //

