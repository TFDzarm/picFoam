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

#include "CountInitialization.H"
#include "polyMesh.H"
#include "Random.H"
#include "meshTools.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CountInitialization<CloudType>::CountInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CountInitialization<CloudType>::~CountInitialization()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::CountInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    const polyMesh& mesh(cloud.mesh());
    Random& rndGen(cloud.rndGen());

    const vector velocity(this->coeffDict().lookup("velocity"));
    const dictionary& countlistDict
    (
        this->coeffDict().subDict("CountList")
    );
    bool countPerCell =  readBool(this->coeffDict().lookup("countPerCell"));

    List<label> countList(cloud.typeIdList().size(),0);
    List<word> molecules(countlistDict.toc());

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);
        label typeId(findIndex(cloud.typeIdList(), moleculeName));
        if(typeId == -1)
            FatalErrorInFunction << "Species " << moleculeName << " is not defined in cloud" << nl << abort(FatalError);

        countList[typeId] = readLabel
        (
            countlistDict.lookup(moleculeName)
        );
    }

    List<label> insertedCount(cloud.typeIdList().size(),0);

    scalar meshVolume = sum(mesh.cellVolumes());
    forAll(mesh.cells(), celli)
    {
        label cellAdd = 0;
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh,
            celli
        );
        scalar cellVolume = mesh.cellVolumes()[celli];


        if(countPerCell) {
            forAll(cloud.typeIdList(), typeId)
            {
                const typename CloudType::parcelType::constantProperties& cP =
                cloud.constProps(typeId);

                label count = countList[typeId];

                for(int i = 0; i < count; i++)
                {
                    List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                    (
                        mesh,
                        celli
                    );
                    label rndTetI = rndGen.sampleAB(0,cellTets.size());
                    const tetIndices& cellTetIs = cellTets[rndTetI];
                    tetPointRef tet = cellTetIs.tet(cloud.mesh());

                    point p = tet.randomPoint(rndGen);
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
        else {
            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];
                tetPointRef tet = cellTetIs.tet(cloud.mesh());
                scalar tetVolume = tet.mag();

                forAll(cloud.typeIdList(), typeId)
                {
                    const typename CloudType::parcelType::constantProperties& cP =
                    cloud.constProps(typeId);

                    scalar count = countList[typeId];

                    // Calculate the number of particles required
                    scalar particlesRequired = count*tetVolume / meshVolume;

                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    label particlesRequiredInCell = label(count*cellVolume / meshVolume);
                    if
                    (
                         (particlesRequired - nParticlesToInsert)
                      > rndGen.scalar01() && cellAdd < particlesRequiredInCell
                    )
                    {
                        nParticlesToInsert++;
                    }

                    insertedCount[typeId] += nParticlesToInsert;
                    cellAdd += nParticlesToInsert;
                    for (label pI = 0; pI < nParticlesToInsert; pI++)
                    {
                        point p = tet.randomPoint(rndGen);
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


    }
    if(countPerCell)
       return;

    forAll(insertedCount, typeId)
    {
        const typename CloudType::parcelType::constantProperties& cP =
        cloud.constProps(typeId);

        label rem = countList[typeId]-insertedCount[typeId];
        for(int i = 0; i < rem; i++)
        {
            label rndCellI = rndGen.sampleAB(0,mesh.cells().size());
            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh,
                rndCellI
            );
            label rndTetI = rndGen.sampleAB(0,cellTets.size());
            const tetIndices& cellTetIs = cellTets[rndTetI];
            tetPointRef tet = cellTetIs.tet(cloud.mesh());

            point p = tet.randomPoint(rndGen);
            meshTools::constrainDirection(cloud.mesh(), cloud.mesh().solutionD(), p);

            vector U = cloud.equipartitionLinearVelocity
            (
                this->temperatures_[typeId],
                cP.mass()
            );

            U += velocity;

            cloud.addNewParcel(p, rndCellI, U, typeId);

        }
    }
}

// ************************************************************************* //

