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

#include "ParticleListInitialization.H"
#include "constants.H"
#include "fvMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleListInitialization<CloudType>::ParticleListInitialization
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InitializationModel<CloudType>(dict, cloud, typeName),
    useElectronVolt_(dict.lookupOrDefault<bool>("electronVolt",false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleListInitialization<CloudType>::~ParticleListInitialization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleListInitialization<CloudType>::initialiseParticles()
{
    CloudType& cloud(this->cloud());
    const fvMesh& mesh(cloud.mesh());

    forAllConstIter(IDLList<entry>, this->coeffDict(), iter)
    {
        if(iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            const vector velocity(subDict.lookup("velocity"));
            label cellI = subDict.lookupOrDefault<label>("cellId",-1);
            vector position;
            if(cellI > -1)
            {
                if(cellI > mesh.C().size())
                    FatalErrorInFunction << "cellId is out if range [0;" << mesh.C().size() << "]" << nl << abort(FatalError);

                position = mesh.C()[cellI];
            }
            else
            {
                position = vector(subDict.lookup("position"));
                cellI = mesh.findCell(position);
                if(cellI == -1)
                    FatalErrorInFunction << "Position outside of mesh" << nl << abort(FatalError);
            }

            const word species(subDict.lookup("species"));

            label typeId = findIndex(cloud.typeIdList(),species);
            if(typeId == -1)
                FatalErrorInFunction << "Unkown species " << species << nl << abort(FatalError);


            const typename CloudType::parcelType::constantProperties& cP =
            cloud.constProps(typeId);

            scalar temperature(readScalar(subDict.lookup("temperature")));
            if(useElectronVolt_)
                temperature *= (constant::electromagnetic::e.value()/constant::physicoChemical::k.value());

            vector pVel = cloud.equipartitionLinearVelocity(temperature,cP.mass());

            pVel += velocity;

            cloud.addNewParcel(position, cellI, pVel, typeId);
        }
    }


}

// ************************************************************************* //

