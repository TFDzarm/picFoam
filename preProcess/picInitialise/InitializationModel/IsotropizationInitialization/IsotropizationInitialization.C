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

#include "IsotropizationInitialization.H"
#include "vector.H"
#include "Random.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropizationInitialization<CloudType>::IsotropizationInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropizationInitialization<CloudType>::~IsotropizationInitialization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IsotropizationInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    Random& rndGen(cloud.rndGen());

    scalar fraction(readScalar(this->coeffDict().lookup("Tperp_over_Tpara")));

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

                    scalar Upara = Foam::sqrt(physicoChemical::k.value()*this->temperatures_[typeId]/(((2.0/3.0)*fraction+(1.0/3.0))*cP.mass()));
                    scalar Uperp = Foam::sqrt(physicoChemical::k.value()*this->temperatures_[typeId]/(((2.0/3.0)+(1.0/(3.0*fraction)))*cP.mass()));



                    vector U = rndGen.sampleNormal<vector>();
                    U[0] *= Upara;
                    U[1] *= Uperp;
                    U[2] *= Uperp;


                    cloud.addNewParcel(p, celli, U, typeId);
                }
            }
        }
    }
}

// ************************************************************************* //

