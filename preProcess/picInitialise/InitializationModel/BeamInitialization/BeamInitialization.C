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

#include "BeamInitialization.H"
#include "meshTools.H"
#include "constants.H"
#include "polyMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BeamInitialization<CloudType>::BeamInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BeamInitialization<CloudType>::~BeamInitialization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::BeamInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    const polyMesh& mesh(cloud.mesh());
    const word speciesA(this->coeffDict().lookup("speciesA"));
    const word speciesB(this->coeffDict().lookup("speciesB"));

    const vector velA(this->coeffDict().lookup("velocityA"));
    const vector velB(this->coeffDict().lookup("velocityB"));

    wordList species;species.append(speciesA);species.append(speciesB);

    const vector mins(this->coeffDict().lookupOrDefault("min",mesh.bounds().min()));
    const vector maxs(this->coeffDict().lookupOrDefault("max",mesh.bounds().max()));
    const vector particleNumbers(this->coeffDict().lookup("particleNumbers"));

    vector displacementDir(this->coeffDict().lookup("displacementDir"));
    const scalarList displacement(this->coeffDict().lookup("displacement"));
    const scalar cycles(readScalar(this->coeffDict().lookup("cycles")));
    //make sure its normalized
    displacementDir /= mag(displacementDir);

    vector span = (maxs-mins);

    scalar spanLength = span&displacementDir;
    vector delta = span;
    delta.x() /= Foam::max(label(particleNumbers.x())-1,2);
    delta.y() /= Foam::max(label(particleNumbers.y())-1,2);
    delta.z() /= Foam::max(label(particleNumbers.z())-1,2);

    forAll(species, spI)
    {
        const word& moleculeName(species[spI]);
        label typeId(findIndex(cloud.typeIdList(), moleculeName));
        if(typeId == -1)
            FatalErrorInFunction << "Species " << moleculeName << " is not defined in cloud" << nl << FatalErrorInFunction;

        for(label i = 0; i < particleNumbers.x();i++)
            for(label j = 0; j < particleNumbers.y(); j++)
                for(label k = 0; k < particleNumbers.z(); k++)
                {
                    vector p(mins.x()+delta.x()*i,mins.y()+delta.y()*j,mins.z()+delta.z()*k);

                    scalar distance = (p-mins)&displacementDir;
                    scalar curDisplacement = displacement[spI]*::sin(distance/spanLength*cycles*constant::mathematical::pi);

                    vector displacedP = p+displacementDir*curDisplacement;

                    const typename CloudType::parcelType::constantProperties& cP =
                    cloud.constProps(typeId);

                    vector U = cloud.equipartitionLinearVelocity
                    (
                        this->temperatures_[typeId],
                        cP.mass()
                    );
                    if(spI == 0)
                        U+= velA;
                    else
                        U+= velB;

                    label cellofPos = mesh.findCell(displacedP);
                    if(cellofPos < 0)
                        FatalErrorInFunction << "Tried to create particle outside of mesh " << displacedP << nl << abort(FatalError);

                    cloud.addNewParcel(displacedP, cellofPos, U, typeId);
                }
    }

}

// ************************************************************************* //

